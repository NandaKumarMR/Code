#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 19:19:46 2020

@author: nandhu
"""

import streamlit as st
#import os
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
smiles = [user_input, '123']
df = pd.DataFrame(smiles)
df.to_csv('moleculeinput.smi', sep='\t', index=False, header=False)
import subprocess

#os.system("cd /home/nandhu/Desktop/datascience/Covid/input/")
subprocess.run(["bash","./padel.sh"])
#process = subprocess.run('padel.sh', shell=True, check=True, timeout=10) 
#import seaborn as sns
#from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
import pickle
dataset = pd.read_csv('descriptors_output.csv')
X = dataset.drop('Name', axis=1)
load_model = pickle.load(open('covid.pkl','rb'))
prediction = load_model.predict(X)
st.write("""
### Predicited pIC50 value is """)
st.write(str(prediction)+"Molar")