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
    except UnboundLocalError:
        multi_mat_1 = pd.DataFrame()

    # Extract Phase 2
    try:
        multi_mat_2 = parse(ds_path_Phase2, rid=row_names_lm,cid = sig_id)
        multi_mat_2 = multi_mat_2.data_df.sort_index()
    except UnboundLocalError:
        multi_mat_2 = pd.DataFrame()


