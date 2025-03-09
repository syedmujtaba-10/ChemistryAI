import os
import time
import requests
import base64
import supabase
from io import BytesIO
from dotenv import load_dotenv
from fastapi import FastAPI, Query, HTTPException, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.security import OAuth2PasswordBearer, OAuth2PasswordRequestForm
from rxn4chemistry import RXN4ChemistryWrapper
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdChemReactions import ChemicalReaction
from langchain_community.chat_models import ChatOpenAI
from langchain_openai import OpenAIEmbeddings
from langchain_community.vectorstores import FAISS
from langchain.chains import RetrievalQA

# Load environment variables
load_dotenv()

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ✅ Load API Keys
IBM_RXN_API_KEY = os.getenv("IBM_RXN_API_KEY")
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_ANON_KEY = os.getenv("SUPABASE_ANON_KEY")

if not IBM_RXN_API_KEY or not OPENAI_API_KEY:
    raise ValueError("⚠️ Please set the IBM_RXN_API_KEY and OPENAI_API_KEY in your environment variables.")

if not SUPABASE_URL or not SUPABASE_ANON_KEY:
    raise ValueError("⚠️ Please set the SUPABASE_URL and SUPABASE_ANON_KEY in your environment variables.")

# ✅ Initialize Supabase Client
supabase_client = supabase.create_client(SUPABASE_URL, SUPABASE_ANON_KEY)

# ✅ Initialize IBM RXN Wrapper
rxn_wrapper = RXN4ChemistryWrapper(api_key=IBM_RXN_API_KEY)
rxn_wrapper.create_project("chemistry_prediction_project")

# ✅ Load FAISS Database
embedding_model = OpenAIEmbeddings(openai_api_key=OPENAI_API_KEY)
vector_store = FAISS.load_local("faiss_atkins_db", embeddings=embedding_model, allow_dangerous_deserialization=True)
retriever = vector_store.as_retriever()

# ✅ Define RAG pipeline using GPT-4
rag_chain = RetrievalQA.from_chain_type(
    llm=ChatOpenAI(model_name="gpt-4", temperature=0.5, openai_api_key=OPENAI_API_KEY),
    retriever=retriever
)

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="auth/login")

# --------------------------------------------------------
# 🔹 AUTHENTICATION ENDPOINTS
# --------------------------------------------------------

@app.post("/auth/register")
async def register_user(email: str, password: str):
    """ Register a new user in Supabase """
    response = supabase_client.auth.sign_up({"email": email, "password": password})
    
    if response.get("error"):
        raise HTTPException(status_code=400, detail=response["error"]["message"])
    
    return {"message": "User registered successfully!", "user": response["user"]}

@app.post("/auth/login")
async def login_user(form_data: OAuth2PasswordRequestForm = Depends()):
    """ Login user and return authentication token """
    response = supabase_client.auth.sign_in_with_password(
        {"email": form_data.username, "password": form_data.password}
    )

    if response.get("error"):
        raise HTTPException(status_code=400, detail=response["error"]["message"])

    return {"access_token": response["session"]["access_token"], "token_type": "bearer"}

@app.get("/auth/user")
async def get_user(token: str = Depends(oauth2_scheme)):
    """ Get user details from token """
    user = supabase_client.auth.get_user(token)
    
    if user.get("error"):
        raise HTTPException(status_code=401, detail="Invalid token")

    return user["user"]

@app.post("/auth/logout")
async def logout_user(token: str = Depends(oauth2_scheme)):
    """ Logout the current user """
    response = supabase_client.auth.sign_out(token)

    if response.get("error"):
        raise HTTPException(status_code=400, detail=response["error"]["message"])

    return {"message": "User logged out successfully"}

# --------------------------------------------------------
# 🔹 CHEMISTRY PREDICTION ENDPOINTS
# --------------------------------------------------------

@app.get("/chemistry/full-explanation")
async def get_full_chemistry_explanation(reactants: str = Query(..., description="Enter reactants (comma-separated names or SMILES)")):
    """
    API Endpoint:
      1. Predicts the chemical product using IBM RXN.
      2. Converts the SMILES reaction to a human-readable equation using GPT-4.
      3. Retrieves a detailed explanation from the RAG agent.
      4. Generates a reaction image using RDKit.
    """
    prediction_result = predict_reaction_product(reactants)
    if "error" in prediction_result:
        return prediction_result

    reaction_smiles = prediction_result["predicted_product_smiles"]
    human_readable_reaction = convert_smiles_to_text(reaction_smiles)
    detailed_explanation = ask_rag_agent(human_readable_reaction)
    reaction_image = generate_reaction_image(reaction_smiles, prediction_result["reactants_smiles"])

    return {
        "reactants": prediction_result["reactants"],
        "reactants_smiles": prediction_result["reactants_smiles"],
        "predicted_product_smiles": prediction_result["predicted_product_smiles"],
        "human_readable_reaction": human_readable_reaction,
        "detailed_explanation": detailed_explanation,
        "reaction_image": reaction_image
    }

@app.get("/chemistry/ask-rag")
async def ask_rag(question: str = Query(..., description="Ask a chemistry-related question")):
    """ Passes a question to the RAG agent for explanation """
    detailed_explanation = ask_rag_agent(question)
    return {"question": question, "detailed_explanation": detailed_explanation}

# --------------------------------------------------------
# 🔹 HELPER FUNCTIONS
# --------------------------------------------------------

def predict_reaction_product(reactants: str):
    """ Predicts reaction product using IBM RXN """
    reactant_list = reactants.split(",")
    smiles_list = [get_smiles(r.strip()) if not r[0].isalpha() else r.strip() for r in reactant_list]
    reactants_smiles = ".".join(smiles_list)
    response = rxn_wrapper.predict_reaction(reactants_smiles)
    
    if "prediction_id" not in response:
        return {"error": "Failed to submit reaction for prediction."}

    prediction_id = response["prediction_id"]
    for _ in range(10):
        time.sleep(2)
        results = rxn_wrapper.get_predict_reaction_results(prediction_id)
        if results["response"]["payload"]["status"] == "SUCCESS":
            return {
                "reactants": reactant_list,
                "reactants_smiles": reactants_smiles,
                "predicted_product_smiles": results["response"]["payload"]["attempts"][0]["smiles"]
            }
    
    return {"error": "IBM RXN prediction took too long or failed."}

def convert_smiles_to_text(smiles_reaction: str):
    """ Converts SMILES reaction to human-readable equation using GPT-4 """
    response = requests.post(
        "https://api.openai.com/v1/chat/completions",
        json={"model": "gpt-4", "messages": [{"role": "user", "content": f"Convert this reaction: {smiles_reaction}"}]},
        headers={"Authorization": f"Bearer {OPENAI_API_KEY}"}
    )
    return response.json()["choices"][0]["message"]["content"].strip() if response.status_code == 200 else "Error converting SMILES."

