# 🔬 **Chemistry AI - Intelligent Chemistry Assistant**  

🚀 **An AI-powered web application for chemistry students and researchers**.  
Predict chemical reactions, ask complex chemistry questions, and visualize reactions – all in one place!  

🧠 **Powered by Retrieval-Augmented Generation (RAG) using Atkins' Physical Chemistry textbook** to provide **accurate, textbook-backed explanations**.  

---

## 🌟 **Features**
### 🔹 **1. Predict Chemical Reactions**
- Input **reactants** (chemical names).
- Get **predicted products** using **IBM RXN for Chemistry API**.
- Convert **SMILES notation** into **human-readable chemical equations** using **GPT-4**.
- **Reaction visualization** with **RDKit**.

### 🔹 **2. AI-Powered Chemistry Q&A (RAG Agent)**
- **Backed by Atkins’ Physical Chemistry textbook**.
- **RAG (Retrieval-Augmented Generation) model** extracts **accurate information** directly from the textbook.
- Provides **detailed, structured, and scientific** explanations for chemistry-related questions.

### 🔹 **3. Interactive Reaction Visualization**
- Displays **reaction diagrams** using **RDKit**.
- Highlights **reactants and products** in color.

### 🔹 **4. User Authentication & Profiles**
- Secure **login & register** with **Supabase**.
- Users can **log in, view their profile, and log out**.

### 🔹 **5. Flashcards for Chemistry Learning**
- **Save** chemistry flashcards with questions and answers.
- **Retrieve** saved flashcards based on **user authentication**.
- **Designed for personalized learning**.

---

## 🎯 **Tech Stack**
- **Frontend**: React + TypeScript + TailwindCSS + Vite  
- **Backend**: FastAPI + Python  
- **Database**: Supabase (PostgreSQL)  
- **AI**:
  - **IBM RXN for Chemistry** (Reaction Prediction)
  - **GPT-4 (OpenAI API)** (SMILES Conversion & Explanations)
  - **RDKit** (Reaction Visualization)
  - **FAISS** (Vector Search for RAG Agent)
  - **RAG Model Trained on Atkins' Physical Chemistry**

---

## 🛠 **Setup & Installation**
### ✅ **1. Clone the Repository**
```sh
git clone https://github.com/your-username/chemistry-ai.git
cd ChemistryAI
```

### ✅ **2. Backend Setup (FastAPI)**
```sh
cd Backend
python -m venv chemenv
source chemenv/bin/activate  # On Windows: chemenv\Scripts\activate
pip install -r requirements.txt
```
**Create a `.env` file** inside `Backend/`:
```
IBM_RXN_API_KEY=your_ibm_rnx_key
OPENAI_API_KEY=your_openai_key
SUPABASE_URL=your_supabase_url
SUPABASE_ANON_KEY=your_supabase_anon_key
```
Then, run the FastAPI server:
```sh
uvicorn main:app --reload
```

---

### ✅ **3. Frontend Setup (React + Vite)**
```sh
cd ../Frontend
npm install
```
**Create a `.env` file** inside `Frontend/`:
```
VITE_SUPABASE_URL=your_supabase_url
VITE_SUPABASE_ANON_KEY=your_supabase_anon_key
VITE_BACKEND_URL=http://127.0.0.1:8000  # Change this if deployed
```
Start the development server:
```sh
npm run dev
```
---


## 🎯 **Usage Guide**
### 🔹 **Predict a Reaction**
1. Enter **reactants** (e.g., "acetic acid, ethanol").
2. Click **"Predict Reaction"**.
3. View:
   - **Human-readable reaction**.
   - **AI explanation**.
   - **Reaction visualization**.

### 🔹 **Ask the RAG Agent (Powered by Atkins' Physical Chemistry)**
1. Enter a **chemistry question**.
2. AI fetches the **most relevant textbook-based explanation**.
3. Get a **detailed, structured response** with **scientific accuracy**.

### 🔹 **Predict Reagents**
1. Enter **starting material and desired product**.
2. Get **reagents list** + **reaction explanation**.

### 🔹 **Flashcards**
1. **Save** important chemistry flashcards.
2. **Retrieve** saved flashcards based on your **user authentication**.

### 🔹 **Login/Register**
- **Sign up/login** using **Supabase**.
- View **your profile** and logout.

---

🔥 **Star the repo if you like this project!** ⭐  
🚀 **Let’s revolutionize chemistry learning with AI!** 🧪
