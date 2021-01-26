import streamlit as st
st.write("""
# Welcome to pIC50 Application

- NANDA KUMAR M R

""")
st.write("""Target: **Coronavirus**""")
st.write("""Organism: **Homo Sapines**""")
st.write("""Target Type: **SINGLE PROTEIN** """)	


st.sidebar.header('User Input Smiles')
user_input = st.sidebar.text_input("Smiles","C=O")
import pandas as pd
mol = str(user_input)

if(not (mol and mol.strip())): 
    st.write("""### Please enter Valid Smiles """) 
else : 
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski

    try:

        desc_MolWt = Descriptors.MolWt(Chem.MolFromSmiles(mol))

        desc_MolLogP = Descriptors.MolLogP(Chem.MolFromSmiles(mol))
        desc_NumHDonors = Lipinski.NumHDonors(Chem.MolFromSmiles(mol))
        desc_NumHAcceptors = Lipinski.NumHAcceptors(Chem.MolFromSmiles(mol))



        values = np.array([desc_MolWt,desc_MolLogP, desc_NumHDonors,desc_NumHAcceptors]).reshape(1,4)
        columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]   

        df_lipinski = pd.DataFrame(data=values,columns=columnNames)

        if df_lipinski.empty:
             st.write("""
            ### Please enter Valid Smiles """)     
        else:
            st.write("""
            ### Lipinski's Descriptors for the Given Smiles """) 
            st.write(df_lipinski)

            from sklearn.ensemble import RandomForestRegressor
            import pickle

            load_model = pickle.load(open('covid_model.pkl','rb'))
            prediction = load_model.predict(df_lipinski)
            st.write("""
            ### Predicited pIC50 value is """)
            st.write(str(prediction)+"Molar")
            if  prediction>=6:
                st.write("""### Predicted State """)
                st.write("Active")
            elif prediction<=5:
                st.write("""### Predicted State """)
                st.write("Inactive")
            else:
                st.write("""### Predicted State """)
                st.write("Intermediate")
            
                
            
    except:
            st.write("""
            ### Please enter Valid Smiles """)

