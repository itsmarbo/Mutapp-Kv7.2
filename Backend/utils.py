#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix
import matplotlib.pyplot as plt
import os

########################################################################################################################

def challenge_encoding(X_train, X_ch):
    
    """
    Description
    -----------
    Do a categorical encoding
    
    Arguments
    ----------
    X_train (numpy array): training data with which we train the model.
    They are needed as a basis in the categorical encoding. 
    Can be the output of sklearn.train_test_split() function.
    
    X_ch (numpy array): data from which we want to predict its label.
        
    Returns
    ----------
    X_ch_enc (numpy array): encoded codified numpy array that models need.
    X_ch_df (DataFrame): encoded codified df. 
    
    """

    # Convert numpy arrays in df
    c = ['affected_domain', 'initial_aa', 'final_aa', 'functional_domain',
         'd_size', 'd_hf', 'd_vol', 'd_msa', 'd_charge', 'd_pol', 'd_aro', 
         'residue_conserv', 'secondary_str', 'plddt', 'str_pos', 'MTR']
    
    X_train_df = pd.DataFrame(X_train, columns = c)
    X_ch_df = pd.DataFrame(X_ch, columns = c)
    
    l = ['affected_domain','initial_aa', 'final_aa', 'functional_domain',
         'd_charge', 'd_pol', 'd_aro', 'secondary_str', 'str_pos']
    
    for i in l: 
        # First: LabelEncoder()
        le = LabelEncoder()
        le.fit(X_train_df[i]) # fit only in training set
        domain_labels_train = le.transform(X_train_df[i])
        domain_labels_ch = le.transform(X_ch_df[i])
        X_train_df["label_train"] = domain_labels_train
        X_ch_df["label_ch"] = domain_labels_ch

        # Second: OneHotEncoder
        oe = OneHotEncoder()
        oe.fit(X_train_df[[i]])
        domain_feature_arr_train = oe.transform(X_train_df[[i]]).toarray()
        domain_feature_arr_ch = oe.transform(X_ch_df[[i]]).toarray()

        domain_feature_labels = list(le.classes_)
        domain_features_train = pd.DataFrame(domain_feature_arr_train, columns = domain_feature_labels)
        domain_features_ch = pd.DataFrame(domain_feature_arr_ch, columns = domain_feature_labels)


        # Update df
        X_ch_df = pd.concat([X_ch_df, domain_features_ch], axis = 1)
        X_ch_df = X_ch_df.drop(["label_ch",i], axis =1)

    return X_ch_df.to_numpy(), X_ch_df
    
    
########################################################################################################################

    
def categorical_encoding(X_train, X_test):
    """
    Description
    -----------
    Do a categorical encoding
    
    Arguments
    ----------
    X_train (numpy array): training data with which we train the model.
    They are needed as a basis in the categorical encoding. 
    Can be the output of sklearn.train_test_split() function.
    
    X_test (numpy array): data from which we want to predict its label.
    Can be the output of sklearn.train_test_split() function.
        
    Returns
    ----------
    X_train_enc (numpy array): encoded codified numpy array that models need.
    X_train_df (DataFrame): encoded codified df. 
    X_test_enc (numpy array): encoded codified numpy array that models need.
    X_test_df (DataFrame): encoded codified df. 
    
    """

    # Convert numpy arrays in df
    c = ['affected_domain', 'initial_aa', 'final_aa', 'functional_domain','d_size', 'd_hf', 'd_vol', 
         'd_msa', 'd_charge', 'd_pol', 'd_aro', 'residue_conserv', 'secondary_str', 'plddt', 'str_pos', 'MTR']
    
    X_train_df = pd.DataFrame(X_train, columns = c)
    X_test_df = pd.DataFrame(X_test, columns = c)
    
    l = ['affected_domain','initial_aa', 'final_aa', 'functional_domain',
         'd_charge', 'd_pol', 'd_aro', 'secondary_str', 'str_pos']
    
    for i in l: 
        # First: LabelEncoder()
        le = LabelEncoder()
        le.fit(X_train_df[i]) # fit only in training set
        domain_labels_train = le.transform(X_train_df[i])
        domain_labels_test = le.transform(X_test_df[i])
        X_train_df["label_train"] = domain_labels_train
        X_test_df["label_test"] = domain_labels_test

        # Second: OneHotEncoder
        oe = OneHotEncoder()
        oe.fit(X_train_df[[i]])
        domain_feature_arr_train = oe.transform(X_train_df[[i]]).toarray()
        domain_feature_arr_test = oe.transform(X_test_df[[i]]).toarray()
        domain_feature_labels = list(le.classes_)
        domain_features_train = pd.DataFrame(domain_feature_arr_train, columns = domain_feature_labels)
        domain_features_test = pd.DataFrame(domain_feature_arr_test, columns = domain_feature_labels)

        # Update df
        X_train_df = pd.concat([X_train_df, domain_features_train], axis = 1)
        X_train_df = X_train_df.drop(["label_train",i], axis =1)
        
        X_test_df = pd.concat([X_test_df, domain_features_test], axis = 1)
        X_test_df = X_test_df.drop(["label_test",i], axis =1)

    
    return X_train_df.to_numpy(), X_test_df.to_numpy(), X_train_df, X_test_df 


########################################################################################################################
  
    


########################################################################################################################

def select_features(X_train, X_train_df, y_train, X_test, n):
    """
  Arguments
  ----------
  X_train: training set.
  y_train: labels of training set.
  X_test: test set.
  n: number of features

  Return
  ---------
  Apply best FS seen for model
    """
    l_feat = []
    l_pos = []
    fs = SelectKBest(score_func=f_classif, k=n)
    fs.fit(X_train, y_train)
    index = fs.fit(X_train, y_train).get_support(indices = True)
    l_pos.append(index)
    l_feat.append(X_train_df.columns[:][index])

    X_train_fs = fs.transform(X_train)
    X_test_fs = fs.transform(X_test)
    
    return X_train_fs, X_test_fs, l_feat, l_pos
    
########################################################################################################################


########################################################################################################################


	
########################################################################################################################
	
def pdb2plddt(in_file):

    # List to store plddt
    plddt = []

    # extract plddt values
    f = open(in_file, "r")
    for l in f:
        if "CA" in l:
            plddt.append(float(l.strip().split()[-2]))
    f.close()

    # return plddt
    return plddt
	
########################################################################################################################

def KCNQ2_DDBB_generation(df):
	
    """
    Arguments
    ----------
    df: single column dataframe with the nomenclature of the mutations for which 
    	descriptors are desired.
    
    Returns
    -----------
    df: dataframe (n rows x n descriptors) with all initial descriptors.
    
    Comments
    -----------
    There are KCNQ2 specific descriptors such as the evolutionary conservation of 
    its residues (residue_conserv), its secondary structure (secondary_str) as well as 
    the functional (functional_domain), topological domains affected ("Location")...that 
    are specific to KCNQ2. If a different channel wants to be used, these descriptors 
    must be modified. 
    
    For the function to work, a file 'FINAL_plot_conservation.csv' with the residue
    conservation values is needed in the same folder where this .py is executed.
    This file will be opened automatically inside the function. 
    """
    
    # initial_aa, final_aa, position_aaa
    df['initial_aa'] = df['Mutationppt'].apply(lambda x: x[0])
    df['final_aa'] = df['Mutationppt'].apply(lambda x: x[-1] if x[-3:] != "del" else "-")
    df['position_aa'] = df["Mutationppt"].apply(lambda x: int(x[1:-1]) if x[-3:] != "del" else int(x[1:-3]))
    df = df.sort_values(by=['position_aa'])
    
    
    # "Location"
    df.loc[(df["position_aa"] >= 1) & (df["position_aa"] <= 91), "Location"] = "Cytoplasmic"
    df.loc[(df["position_aa"] >= 92) & (df["position_aa"] <= 112), "Location"] = "S1"
    df.loc[(df["position_aa"] >= 113) & (df["position_aa"] <= 122), "Location"] = "Extracelullar"
    df.loc[(df["position_aa"] >= 123) & (df["position_aa"] <= 143), "Location"] = "S2"
    df.loc[(df["position_aa"] >= 144) & (df["position_aa"] <= 166), "Location"] = "Cytoplasmic"
    df.loc[(df["position_aa"] >= 167) & (df["position_aa"] <= 187), "Location"] = "S3"
    df.loc[(df["position_aa"] >= 188) & (df["position_aa"] <= 195), "Location"] = "Extracelullar"
    df.loc[(df["position_aa"] >= 196) & (df["position_aa"] <= 218), "Location"] = "S4"
    df.loc[(df["position_aa"] >= 219) & (df["position_aa"] <= 231), "Location"] = "Cytoplasmic"
    df.loc[(df["position_aa"] >= 232) & (df["position_aa"] <= 252), "Location"] = "S5"
    df.loc[(df["position_aa"] >= 253) & (df["position_aa"] <= 264), "Location"] = "Extracelullar"
    df.loc[(df["position_aa"] >= 265) & (df["position_aa"] <= 285), "Location"] = "Pore"
    df.loc[(df["position_aa"] >= 286) & (df["position_aa"] <= 291), "Location"] = "Extracelullar"
    df.loc[(df["position_aa"] >= 292) & (df["position_aa"] <= 312), "Location"] = "S6"
    df.loc[(df["position_aa"] >= 313) & (df["position_aa"] <= 331), "Location"] = "Cytoplasmic"
    df.loc[(df["position_aa"] >= 332) & (df["position_aa"] <= 350), "Location"] = "hA"
    df.loc[(df["position_aa"] >= 351) & (df["position_aa"] <= 356), "Location"] = "Cytoplasmic"
    df.loc[(df["position_aa"] >= 357) & (df["position_aa"] <= 366), "Location"] = "hTW"
    df.loc[(df["position_aa"] >= 367) & (df["position_aa"] <= 534), "Location"] = "Cytoplasmic"
    df.loc[(df["position_aa"] >= 535) & (df["position_aa"] <= 559), "Location"] = "hB"
    df.loc[(df["position_aa"] >= 560) & (df["position_aa"] <= 562), "Location"] = "Cytoplasmic"
    df.loc[(df["position_aa"] >= 563) & (df["position_aa"] <= 594), "Location"] = "hC"
    df.loc[(df["position_aa"] >= 595) & (df["position_aa"] <= 621), "Location"] = "Cytoplasmic"
    df.loc[(df["position_aa"] >= 622) & (df["position_aa"] <= 647), "Location"] = "hD"
    df.loc[(df["position_aa"] >= 648) & (df["position_aa"] <= 872), "Location"] = "Cytoplasmic"
    
    
    # "functional_domain"
    df.loc[(df["position_aa"] >= 1) & (df["position_aa"] <= 91), "functional_domain"] = "unknown_function" 
    df.loc[(df["position_aa"] >= 92) & (df["position_aa"] <= 218), "functional_domain"] = "voltage_domain" 
    df.loc[(df["position_aa"] >= 219) & (df["position_aa"] <= 231), "functional_domain"] = "unknown_function" 
    df.loc[(df["position_aa"] >= 232) & (df["position_aa"] <= 276), "functional_domain"] = "pore_domain" 
    df.loc[(df["position_aa"] >= 277) & (df["position_aa"] <= 282), "functional_domain"] = "selectivity_filter" 
    df.loc[(df["position_aa"] >= 283) & (df["position_aa"] <= 312), "functional_domain"] = "pore_domain" 
    df.loc[(df["position_aa"] >= 313) & (df["position_aa"] <= 331), "functional_domain"] = "unknown_function" 
    df.loc[(df["position_aa"] >= 332) & (df["position_aa"] <= 559), "functional_domain"] = "CaM_interaction" 
    df.loc[(df["position_aa"] >= 560) & (df["position_aa"] <= 562), "functional_domain"] = "unknown_function" 
    df.loc[(df["position_aa"] >= 563) & (df["position_aa"] <= 647), "functional_domain"] = "SID_domain" 
    df.loc[(df["position_aa"] >= 648) & (df["position_aa"] <= 872), "functional_domain"] = "unknown_function"
    
    
    # char_initial_aa, char_final_aa
    df['char_initial_aa'] = df['initial_aa']
    df['char_final_aa'] = df['final_aa']
    
    charge_map = { 'S':'neutral', 'T':'neutral', 'Q':'neutral', 'N':'neutral', 'Y':'neutral', 'C':'neutral', 'G':'neutral', 'A':'neutral', 'V':'neutral', 'L':'neutral', 'I':'neutral', 'M':'neutral', 'P':'neutral', 'F':'neutral', 'W':'neutral', 'D':'negative_acidic', 'E':'negative_acidic', 'K':'positive_basic', 'R':'positive_basic', 'H':'positive_basic' }
    
    df["char_initial_aa"].replace(charge_map, inplace=True)
    df["char_final_aa"].replace(charge_map, inplace=True) 
    
    
    # pol_initial_aa, pol_final_aa
    df['pol_initial_aa'] = df['initial_aa']
    df['pol_final_aa'] = df['final_aa']
    
    pol_map = { 'S':'polar', 'T':'polar', 'Q':'polar', 'N':'polar', 'Y':'non_polar', 'C':'polar', 'G':'polar', 'A':'non_polar', 'V':'non_polar', 'L':'non_polar', 'I':'non_polar', 'M':'non_polar', 'P':'non_polar', 'F':'non_polar', 'W':'non_polar', 'D':'polar', 'E':'polar', 'K':'polar', 'R':'polar', 'H':'polar', }
    
    df["pol_initial_aa"].replace(pol_map, inplace=True) 
    df["pol_final_aa"].replace(pol_map, inplace=True)
    
    
    # aro_final_aa, aro_initial_aa
    df['aro_initial_aa'] = df['initial_aa']
    df['aro_final_aa'] = df['final_aa']
    
    aro_map = { 'S':'non_aromatic', 'T':'non_aromatic', 'Q':'non_aromatic', 'N':'non_aromatic', 'Y':'aromatic', 'C':'non_aromatic', 'G':'non_aromatic', 'A':'non_aromatic', 'V':'non_aromatic', 'L':'non_aromatic', 'I':'non_aromatic', 'M':'non_aromatic', 'P':'non_aromatic', 'F':'aromatic', 'W':'aromatic', 'D':'non_aromatic', 'E':'non_aromatic', 'K':'non_aromatic', 'R':'non_aromatic', 'H':'non_aromatic' }
    
    df["aro_initial_aa"].replace(aro_map, inplace=True)
    df["aro_final_aa"].replace(aro_map, inplace=True)
    
    
    #mw_initial_aa, mw_final_aa 
    df['mw_initial_aa'] = df['initial_aa'] 
    df['mw_final_aa'] = df['final_aa'] 
    
    mw_map = { 'S': 105.09, 'T': 119.12, 'Q': 146.15, 'N': 132.12, 'Y': 181.19, 'C': 121.16, 'G': 75.07, 'A': 89.09, 'V': 117.15, 'L': 131.17, 'I': 131.17, 'M': 149.21, 'P': 115.13, 'F': 165.19, 'W': 204.23, 'D': 133.10 , 'E': 147.13, 'K': 146.19, 'R': 174.20, 'H': 155.16, } 
    
    df["mw_initial_aa"].replace(mw_map, inplace=True) 
    df["mw_final_aa"].replace(mw_map, inplace=True)
    
    
    # v_e_initial_aa, v_e_final_a
    df['v_e_initial_aa'] = df['initial_aa']
    df['v_e_final_aa'] = df['final_aa']
    
    v_map = { 'K': 68.0, 'H': 49.2, 'R': 70.8, 'D': 31.3, 'E': 47.2, 'N': 35.4, 'Q': 51.3, 'S': 18.1, 'T': 34.0, 'C': 28.0, 'G': 0.0, 'A': 15.9, 'P': 41.0, 'V': 47.7, 'M': 62.8, 'I': 63.6, 'L': 63.6, 'Y': 78.5, 'F': 77.2, 'W': 100.0 }
    
    df["v_e_initial_aa"].replace(v_map, inplace=True)
    df["v_e_final_aa"].replace(v_map, inplace=True)
    
    
    # pol_e_inital_aa, pol_e_final_aa
    df['pol_e_initial_aa'] = df['initial_aa']
    df['pol_e_final_aa'] = df['final_aa']
    
    p_map = { 'K': 64.2, 'H': 43.2, 'R': 51.9, 'D': 100.0, 'E': 93.8, 'N': 63.0, 'Q': 45.7, 'S': 32.1, 'T': 21.0, 'C': 7.4, 'G': 37.0, 'A': 25.9, 'P': 21.0, 'V': 8.6, 'M': 4.9, 'I': 0.0, 'L': 0.0, 'Y': 9.9, 'F': 1.2, 'W': 4.9 }
    
    df["pol_e_initial_aa"].replace(p_map, inplace=True)
    df["pol_e_final_aa"].replace(p_map, inplace=True)
    
    
    #ip_e_initial_aa, ip_e_final_aa
    df['ip_e_initial_aa'] = df['initial_aa']
    df['ip_e_final_aa'] = df['final_aa']
    
    ip_e_map = { 'K': 86.9 , 'H': 59.2, 'R': 100.0, 'D': 0.0, 'E': 3.2, 'N': 31.3, 'Q': 34.4, 'S': 34.8, 'T': 45.7, 'C': 26.3, 'G': 38.5, 'A': 39.2, 'P': 40.2, 'V': 38.5, 'M': 35.7, 'I': 39.2, 'L': 38.6, 'Y': 34.4, 'F': 38.6, 'W': 37.7 }
    
    df["ip_e_initial_aa"].replace(ip_e_map, inplace=True)
    df["ip_e_final_aa"].replace(ip_e_map, inplace=True)
    
    
    # hf_e_initial_aa, hf_e_final_aa
    df['hf_e_initial_aa'] = df['initial_aa']
    df['hf_e_final_aa'] = df['final_aa']
    
    hf_e_map = { 'K': 43.5, 'H': 23.1, 'R': 22.6, 'D': 17.5, 'E': 17.8, 'N': 2.4, 'Q': 0.0, 'S': 1.9, 'T': 1.9, 'C': 40.3, 'G': 2.7, 'A': 23.1, 'P': 73.5, 'V': 49.6, 'M': 44.3, 'I': 83.6, 'L': 57.6, 'Y': 70.8, 'F': 76.1, 'W': 100.0 }
    
    df["hf_e_initial_aa"].replace(hf_e_map, inplace=True)
    df["hf_e_final_aa"].replace(hf_e_map, inplace=True)	
    
    #msa_initial_aa, msa_final_aa
    df['msa_e_initial_aa'] = df['initial_aa']
    df['msa_e_final_aa'] = df['final_aa']
    
    msa_map = { 'K': 54.3, 'H': 28.1, 'R': 50.1, 'D': 45.0, 'E': 48.6, 'N': 46.1, 'Q': 43.6, 'S': 40.5, 'T': 35.3, 'C': 7.4, 'G': 54.0, 'A': 37.4, 'P': 66.2, 'V': 19.6, 'M': 3.9, 'I': 7.5, 'L': 10.1, 'Y': 30.1, 'F': 5.5, 'W': 13.8 }
    
    df["msa_e_initial_aa"].replace(msa_map, inplace=True)
    df["msa_e_final_aa"].replace(msa_map, inplace=True)    
    
    #hf_initial_aa, hf_final_aa
    df['hf_initial_aa'] = df['initial_aa']
    df['hf_final_aa'] = df['final_aa']
    
    hf_map = { 'K': -3.9, 'H': -3.2, 'R': -4.5, 'D': -3.5, 'E': -3.5, 'N': -3.5, 'Q': -3.5, 'S': -0.8, 'T': -0.7, 'C': 2.5, 'G': -0.4, 'A': 1.8, 'P': -1.6, 'V': 4.2, 'M': 1.9, 'I': 4.5, 'L': 3.8, 'Y': -1.3, 'F': 2.8, 'W': -0.9 }
    
    df["hf_initial_aa"].replace(hf_map, inplace=True)
    df["hf_final_aa"].replace(hf_map, inplace=True)
    
    
    # d_size, d_hf, d_vol, d_msa
    df["d_size"] = df["mw_initial_aa"] - df["mw_final_aa"]
    df["d_hf"] = df["hf_initial_aa"] - df["hf_final_aa"]
    df["d_vol"] = df["v_e_initial_aa"] - df["v_e_final_aa"]
    df["d_msa"] = df["msa_e_initial_aa"] - df["msa_e_final_aa"]
    
    # d_charge: 
    df.loc[(df["char_initial_aa"] == "positive_basic") & (df["char_final_aa"] == "positive_basic"), "d_charge"] = "pos_to_pos"
    df.loc[(df["char_initial_aa"] == "positive_basic") & (df["char_final_aa"] == "negative_acidic"), "d_charge"] = "pos_to_neg" 
    df.loc[(df["char_initial_aa"] == "positive_basic") & (df["char_final_aa"] == "neutral"), "d_charge"] = "pos_to_neu" 
    df.loc[(df["char_initial_aa"] == "negative_acidic") & (df["char_final_aa"] == "positive_basic"), "d_charge"] = "neg_to_pos" 
    df.loc[(df["char_initial_aa"] == "negative_acidic") & (df["char_final_aa"] == "negative_acidic"), "d_charge"] = "neg_to_neg" 
    df.loc[(df["char_initial_aa"] == "negative_acidic") & (df["char_final_aa"] == "neutral"), "d_charge"] = "neg_to_neu" 
    df.loc[(df["char_initial_aa"] == "neutral") & (df["char_final_aa"] == "positive_basic"), "d_charge"] = "neu_to_pos" 
    df.loc[(df["char_initial_aa"] == "neutral") & (df["char_final_aa"] == "negative_acidic"), "d_charge"] = "neu_to_neg" 
    df.loc[(df["char_initial_aa"] == "neutral") & (df["char_final_aa"] == "neutral"), "d_charge"] = "neu_to_neu" 
    
    # d_pol: 
    df.loc[(df["pol_initial_aa"] == "polar") & (df["pol_final_aa"] == "polar"), "d_pol"] = "p_to_p" 
    df.loc[(df["pol_initial_aa"] == "polar") & (df["pol_final_aa"] == "non_polar"), "d_pol"] = "p_to_np" 
    df.loc[(df["pol_initial_aa"] == "non_polar") & (df["pol_final_aa"] == "polar"), "d_pol"] = "np_to_p" 
    df.loc[(df["pol_initial_aa"] == "non_polar") & (df["pol_final_aa"] == "non_polar"), "d_pol"] = "np_to_np" 
    
    # d_aro: 
    df.loc[(df["aro_initial_aa"] == "aromatic") & (df["aro_final_aa"] == "aromatic"), "d_aro"] = "a_to_a" 
    df.loc[(df["aro_initial_aa"] == "aromatic") & (df["aro_final_aa"] == "non_aromatic"), "d_aro"] = "a_to_na"
    df.loc[(df["aro_initial_aa"] == "non_aromatic") & (df["aro_final_aa"] == "aromatic"), "d_aro"] = "na_to_a"
    df.loc[(df["aro_initial_aa"] == "non_aromatic") & (df["aro_final_aa"] == "non_aromatic"), "d_aro"] = "na_to_na"
    
    # residue_conserv
    condf = pd.read_csv('FINAL_plot_conservation.csv', header = "infer") 
    d = {} 
    for v in condf.iterrows():
    	res = int(v[1][0]) # residues 
    	con = float(v[1][1]) # conservation 
    	d[res] = con 
    res = [i for i in d.keys()] 
    con = [i for i in d.values()] 
    lpos = [i for i in df["position_aa"]] 
    lcons = [] 
    
    for i in range(len(lpos)): 
    	for j in range(len(res)): 
    		if lpos[i] == res[j]: 
    			lcons.append(con[j]) 
    df["r_conserv"] = lcons
    
    # secondary_str 
    prot2 = "CCCCCCCCCCCCCCHHHHHHHHCCECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCHHHHHHHHHHHHHHCCCCHHHTTTTTTTTTTTTTTTTTTTTTTCHHHHHHHHHCTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCHHHHHHHCTTTTTTTTTTTTTTTTTTCHHHHHHHHHHHHCTTTTTTTTTTTTTTTTCCHHHHHHHHHHHHCCCTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCTTTTTTTTTTTTTCCCCCCCCCCCTTTTTTTTTTTTTTTTTTTCTTTTTTTTCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHCCCCCCCHHHCCCCCCCCCCCCCCCCCCCCCCCCCCCCHHHHHCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCECHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHCHHHHHHHHHHHHHHHHCCCCCCCCCCCCCCCCCCCCCCCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHCCCCCCCCCECCCCCCCCCCCCCCCCCCCCCCCCCCCEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEEECCCCHHHHHHHHHHHCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
   
    lss = []

    for i in range(len(prot2)): 
    	for j in lpos:
    		if i+1 == j:
    			lss.append(prot2[i])
    df["sec_str"] = lss
    
    ss_map = { 'C':'coil', 'T':'membrane_helix', 'H':'helix', 'E':'beta_strand' }
    df["sec_str"].replace(ss_map, inplace=True)
    
    
    # plddt of KCNQ2 AF prediction
    data = [f for f in os.listdir() if f.endswith(".pdb")]
    
    for f in data:
    	plddt = pdb2plddt(f)
    	lpos = [i for i in df["position_aa"]]
    	l_plddt = [] 
    	
    
    for i in range(len(lpos)): 
    	for j in range(len(plddt)): 
    		if lpos[i] == j+1:
    			l_plddt.append(plddt[j]) # Create that column in df:
    			
    df["pLDDT"] = l_plddt
    
    
    # str_pos 
    df.loc[(df["position_aa"] >= 1) & (df["position_aa"] <= 69), "str_pos"] = "GAP1" 
    df.loc[(df["position_aa"] >= 70) & (df["position_aa"] <= 184), "str_pos"] = "FRAG1"
    df.loc[(df["position_aa"] >= 185) & (df["position_aa"] <= 194), "str_pos"] = "GAP2"
    df.loc[(df["position_aa"] >= 195) & (df["position_aa"] <= 367), "str_pos"] = "FRAG2"
    df.loc[(df["position_aa"] >= 368) & (df["position_aa"] <= 534), "str_pos"] = "GAP3"
    df.loc[(df["position_aa"] >= 535) & (df["position_aa"] <= 595), "str_pos"] = "FRAG3"
    df.loc[(df["position_aa"] >= 596) & (df["position_aa"] <= 872), "str_pos"] = "GAP4"
    
    # MTR
    mtr = pd.read_csv('MTR.csv', header = "infer", sep = "\t")
    mtr_l = [i for i in mtr["mtr"]]
    lpos = [i for i in df["position_aa"]]
    lmtr = [] 
    
    for i in lpos: 
    	for j in range(len(mtr_l)): 
    		if i == (j+1):
    			lmtr.append(mtr_l[j])
    df["MTR"] = lmtr
    
    
    return df
    	
########################################################################################################################
def preprocessing_ch(df):
    """
    Description: 
        Takes as input the database with the designed columns. 
        It eliminates those columns that are useless, renames the
        remaining ones and performs a binary coding of the label.
    
    Input: 
        Original dataframe
        
    Return:
        Modified dataframe
    """
    
    # Change column names
    df.columns  = ['mutation', 'initial_aa', 'final_aa', 'position_aa', 'topological_domain',
       'functional_domain', 'char_initial_aa', 'char_final_aa', 'pol_initial_aa', 'pol_final_aa',
       'aro_initial_aa', 'aro_final_aa', 'mw_initial_aa', 'mw_final_aa',
       'v_e_initial_aa', 'v_e_final_aa', 'pol_e_initial_aa', 'pol_e_final_aa',
       'ip_e_initial_aa', 'ip_e_final_aa', 'hf_e_initial_aa', 'hf_e_final_aa',
       'msa_e_initial_aa', 'msa_e_final_aa', 'hf_initial_aa', 'hf_final_aa',
       'd_size', 'd_hf', 'd_vol', 'd_msa','d_charge', 'd_pol', 'd_aro', 'residue_conserv', 'secondary_str', 
        'pLDDT', 'str_pos', 'MTR']
    
    # Delete extra columns 
    df = df.drop(['mutation', 'position_aa','char_initial_aa', 'char_final_aa', 'pol_initial_aa', 'pol_final_aa',
       'aro_initial_aa', 'aro_final_aa', 'mw_initial_aa', 'mw_final_aa',
       'v_e_initial_aa', 'v_e_final_aa', 'pol_e_initial_aa', 'pol_e_final_aa',
       'ip_e_initial_aa', 'ip_e_final_aa', 'hf_e_initial_aa', 'hf_e_final_aa',
       'msa_e_initial_aa', 'msa_e_final_aa', 'hf_initial_aa', 'hf_final_aa'], axis = 1)
    
    return df

