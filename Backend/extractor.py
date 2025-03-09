import os
import faiss
import chromadb
import openai
import numpy as np
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.embeddings import OpenAIEmbeddings  #
from langchain_community.vectorstores import FAISS  #
from pypdf import PdfReader


OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")

PDF_PATH = r"D:\ChemistryRAG\Backend\data\AtkinsPhysicalChemistry.pdf"

def extract_text_from_pdf(pdf_path):
    """Extracts text from a PDF file."""
    reader = PdfReader(pdf_path)
    text = ""
    for page in reader.pages:
        text += page.extract_text() + "\n"
    return text


raw_text = extract_text_from_pdf(PDF_PATH)
text_splitter = RecursiveCharacterTextSplitter(chunk_size=500, chunk_overlap=50)
chunks = text_splitter.split_text(raw_text)

embedding_model = OpenAIEmbeddings(openai_api_key=os.getenv("OPENAI_API_KEY"))


vector_store = FAISS.from_texts(chunks, embedding_model)

vector_store.save_local("faiss_atkins_db")
