#!/usr/bin/env python3

from HRProfiler.scripts import HRProfiler as HR
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import pandas as pd

input_dir = "SigProfilerAssignment/SBS/input/"
project = "HRProfiler"

#manually generate the input matrix for HRProfiler
matrices = matGen.SigProfilerMatrixGeneratorFunc(
	project, "GRCh38", input_dir,plot=True, exome=False, bed_file=None, 
	chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100
)

#SBS matrix
matrices["96"].to_csv("HRProfiler/matrixSBS.txt", sep="\t", index=True)
sbs = pd.read_table("HRProfiler/matrixSBS.txt", sep="\t")
sbs.columns = [col.split("_")[0] for col in sbs.columns]

mask = sbs["MutationType"].str.contains(r"\[C>T\]G")
summed1 = sbs.loc[mask, sbs.columns != "MutationType"].sum()

mask2 = sbs["MutationType"].str.contains(r"\[C>G\]T")
summed2 = sbs.loc[mask2, sbs.columns != "MutationType"].sum()
sbs_df = pd.DataFrame(
    [summed1, summed2],
    index=["NCTG", "NCGT"]  # row labels
).reset_index().rename(columns={"index": "mutation_type"})


#CN matrix
cnv = pd.read_table("SigProfilerAssignment/CN/input/FACETS.CNV48.matrix.tsv", sep="\t")
mask = cnv["MutationType"].str.contains("LOH:1Mb-10Mb|LOH:10Mb-40Mb" )
summed1 = cnv.loc[mask, cnv.columns != "MutationType"].sum()

mask2 = cnv["MutationType"].str.contains("2:het:>40Mb|3-4:het:>40Mb")
summed2 = cnv.loc[mask2, cnv.columns != "MutationType"].sum()

mask3 = cnv["MutationType"].str.contains("3-4:het:10Mb-40Mb|5-8:het:10Mb-40Mb")
summed3 = cnv.loc[mask3, cnv.columns != "MutationType"].sum()

cn_df = pd.DataFrame(
    [summed1, summed2,summed3],
    index=["LOH.1.40Mb", "2-4:HET.40Mb", "3-9:HET.10.40Mb"]  # row labels
).reset_index().rename(columns={"index": "mutation_type"})

#ID matrix
matrices["ID"].to_csv("HRProfiler/matrixID.txt", sep="\t", index=True)
id = pd.read_table("HRProfiler/matrixID.txt", sep="\t")
id.columns = [col.split("_")[0] for col in id.columns]
mask = id["Unnamed: 0"].str.contains("5:Del:M" )
summed1 = id.loc[mask, id.columns != "Unnamed: 0"].sum()
# Combine into new DataFrame
id_df = pd.DataFrame(
    [summed1],
    index=["DEL_5_MH"]  # row labels
).reset_index().rename(columns={"index": "mutation_type"})

#generate final input matrix for HRProfiler
combined = pd.concat([sbs_df, cn_df, id_df], axis=0, ignore_index=True)
df_t = combined.T 
df_t = df_t.reset_index().rename(columns={"index": "sample"})
new_columns = df_t.iloc[0]  
df = df_t[1:]              
df.columns = new_columns
df.columns.values[0] = "samples"
df.to_csv("HRProfiler/data_matrix.txt", sep="\t", index=True)

#run HRProfiler
HR.HRProfiler(data_matrix=df,
              genome='GRCh38', 
              exome=False, 
              INDELS_DIR=None,
              SNV_DIR=None,
              CNV_DIR=None, 
			  RESULT_DIR='HRProfiler/output/',              
			  cnv_file_type='ASCAT',
              bootstrap=False, 
              nreplicates=20,
              normalize=True, 
              hrd_prob_thresh=0.5,
              plot_predictions=True,
              organ='BREAST')