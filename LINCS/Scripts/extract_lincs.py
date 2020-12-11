# coding: utf-8

#Import packages
from __future__ import print_function
import sys
import os
import pandas as pd
import numpy as np
import cmapPy
import h5py
from cmapPy.pandasGEXpress.parse import parse
from pandas.util.testing import assert_frame_equal
import time
from numpy import unravel_index,fill_diagonal,nanargmax,nanargmin,nanmax,nanmin
import random
import scipy
from scipy import stats
from scipy.stats import spearmanr, pearsonr
import math
from math import sqrt
from joblib import Parallel, delayed
import multiprocessing
from datetime import datetime
from optparse import OptionParser
startTime = datetime.now()

# Specify directories for metadata files
data_dir = "../../../LINCS-Extraction/data/"

ds_path_Phase1=data_dir+"GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
ds_path_Phase2=data_dir+"GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
gene_info=pd.read_table(data_dir+"GSE92742_Broad_LINCS_gene_info.txt", sep='\t')
sig_info_1=pd.read_table(data_dir+"GSE92742_Broad_LINCS_sig_info.txt", sep="\t")
sig_info_2=pd.read_table(data_dir+"GSE70138_Broad_LINCS_sig_info.txt", sep="\t")

# Load list of sig_ids
all_sig_ids = pd.read_csv("../Output/sig_ids.csv",header=None)
all_sig_ids = all_sig_ids[0].tolist()
#all_sig_ids = all_sig_ids[0:5] #test

# Initialise empty dataframe
empty_df = pd.DataFrame()

# Fill with gene IDs and sig ids
row_names_lm = list((gene_info['pr_gene_id'][gene_info["pr_is_lm"] == 1]).sort_values())
column_names_lm = list(all_sig_ids)
matrix_test_lm = np.empty((len(row_names_lm),len(all_sig_ids)))
matrix_test_lm[:] = np.NaN
row_names_lm=[str(i) for i in row_names_lm]
df_lm = pd.DataFrame(matrix_test_lm, columns=column_names_lm, index=row_names_lm)
compound_geneExpression_lm=df_lm.sort_index()
headers = row_names_lm
headers = ['Sig_id'] + headers

# Write empty matrix
mat = open('../Output/extracted_signatures.txt','w+')
mat.write('\t'.join(map(str,headers)) + '\n')

# Define function for extracting data
def extractLincs(sig_id):

    # Progress
    print(str(all_sig_ids.index(sig_id)), 'out of', str((len(all_sig_ids))))

    # Extract Phase 1
    try:
        multi_mat_1 = parse(ds_path_Phase1, rid = row_names_lm, cid = sig_id)
        multi_mat_1 = multi_mat_1.data_df.sort_index()
    except:
        multi_mat_1 = pd.DataFrame()

    # Extract Phase 2
    try:
        multi_mat_2 = parse(ds_path_Phase2, rid=row_names_lm,cid = sig_id)
        multi_mat_2 = multi_mat_2.data_df.sort_index()
    except:
        multi_mat_2 = pd.DataFrame()

    if not multi_mat_2.equals(empty_df) and not multi_mat_1.equals(empty_df):
        multi_mat=pd.concat([multi_mat_1, multi_mat_2], axis=1)
    elif not multi_mat_1.equals(empty_df) and multi_mat_2.equals(empty_df):
        multi_mat=multi_mat_1
    elif not multi_mat_2.equals(empty_df) and multi_mat_1.equals(empty_df):
        multi_mat=multi_mat_2
    elif multi_mat_2.equals(empty_df) and multi_mat_1.equals(empty_df):
        multi_mat=pd.DataFrame()

    if multi_mat.shape[1]==1:
        #Save data
        single_sig = multi_mat.median(axis=1)
        single_sig = single_sig.tolist()
        single_sig.insert(0,sig_id)
        single_sig = list(map(str, single_sig))
        return single_sig

    if multi_mat.shape[1]==2:
        #Consensus signature is just median
        consensus_sig = multi_mat.median(axis=1)
        consensus_sig = consensus_sig.tolist()
        consensus_sig.insert(0,sig_id)
        consensus_sig = list(map(str, consensus_sig))
        return consensus_sig

    if multi_mat.shape[1]>2:
        #Get consensus signature
        corr_reps = pd.DataFrame(multi_mat.corr(method='spearman'))
        a = np.array(corr_reps)
        np.fill_diagonal(a,'NaN')
        corr_reps_nan = pd.DataFrame(a)
        col_weights = []
        for col in corr_reps_nan:
            col_weight = corr_reps_nan[col].sum()
            col_weights.append(col_weight)
        col_weights_norm = [float(i)/sum(col_weights) for i in col_weights]
        linear_list = []
        for idx, weight in enumerate(col_weights_norm):
            linear = multi_mat.ix[:,idx]*col_weights_norm[idx]
            linear_list.append(linear)
        combined_sig = sum(linear_list)
        combined_sig = list(map(str, combined_sig))
        consensus_sig = combined_sig
        consensus_sig.insert(0,sig_id)
        return consensus_sig

results = Parallel(n_jobs=20,backend="multiprocessing")(delayed(extractLincs)(sig_id) for sig_id in all_sig_ids)

for sig_id in results:
    mat.write('\t'.join(map(str,sig_id)) + '\n')

mat.close()
print("Finished.")
