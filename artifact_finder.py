import numpy as np
import pandas as pd
from tqdm import tqdm
import pickle
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import pyranges as pr


# columns to load from maf
cols = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode']
# load maf
maf = pd.read_csv('data/data_mutations_extended_7.6-consortium.txt', sep='\t', usecols=cols, low_memory=False)

# load sample table
sample = pd.read_csv('data/data_clinical_sample_7.6-consortium.txt', sep='\t', skiprows=4)

# load bed file
gi = pd.read_csv('data/genomic_information_7.6-consortium.txt', sep='\t', low_memory=False)
beds = dict(list(gi[['Chromosome', 'Start_Position', 'End_Position', 'SEQ_ASSAY_ID']].groupby('SEQ_ASSAY_ID')))

# add assay id to amf
maf = pd.merge(maf, sample[['SAMPLE_ID', 'SEQ_ASSAY_ID']], how='left', left_on='Tumor_Sample_Barcode', right_on='SAMPLE_ID')
# make unique variant string
maf['var_str'] = maf['Chromosome'].astype(str) + '_' + maf['Start_Position'].astype(str) + '_' + maf['End_Position'].astype(str) + '_' + maf['Reference_Allele'].astype(str) + '_' + maf['Tumor_Seq_Allele2'].astype(str)

# counts of mutations across assay ID (nan = not in data)
call_counts = maf.groupby(['var_str', 'SEQ_ASSAY_ID']).size().to_frame().reset_index().pivot(index='var_str', columns='SEQ_ASSAY_ID', values=0)

# get counts of sample per assay
sample_counts = sample.groupby('SEQ_ASSAY_ID').size()

# get unique mutations
var_uniq = maf[['var_str', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']].drop_duplicates().set_index('var_str')
# make PyRange for unique variants
var_pr = pr.PyRanges(chromosomes=var_uniq['Chromosome'], starts=var_uniq['Start_Position'], ends=var_uniq['End_Position'])
# include the variant unique identifier, which is used to index the df
var_pr.var_str = var_uniq.index.values

# find overlap of variants in various bed files and add to unique variant table
for bed_name in beds.keys():
    var_uniq[bed_name] = False
    var_uniq.loc[var_pr.overlap(pr.PyRanges(chromosomes=beds[bed_name]['Chromosome'], starts=beds[bed_name]['Start_Position'], ends=beds[bed_name]['End_Position'])).var_str, bed_name] = True

# is nan to true zeros based on coverage of panel for specific variant
call_counts.values[(call_counts.isna() & var_uniq.iloc[:, 6:]).values] = 0

bed_review_idx = dict()
# criteria for flagging variants for assay bias
for bed_name in beds.keys():
    bed_name_other = np.setdiff1d(list(beds.keys()), bed_name)
    bed_review_idx[bed_name] = np.where((call_counts[bed_name] >= 10) & (call_counts[bed_name_other].sum(axis=1, skipna=True) < 10) & ((~call_counts[bed_name_other].isna()).sum(axis=1) >= 3))[0]

bed_name = 'MSK-IMPACT468'
bed_name = 'JHU-50GP'
bed_name = 'UHN-555-V1'
bed_name = 'PHS-FOCUS-V1'

bed_name_other = np.setdiff1d(list(beds.keys()), bed_name)
p = call_counts.iloc[bed_review_idx[bed_name]].apply(lambda x: pd.Series(fisher_exact([[x[bed_name], sample_counts[bed_name]], [x[bed_name_other].sum(skipna=True), sample_counts[bed_name_other[(~x[bed_name_other].isna())]].sum()]], alternative='greater')[1], index=['p-value']), axis=1)
review_table = pd.concat([var_uniq.loc[call_counts.index[bed_review_idx[bed_name]]].iloc[:, :6], call_counts[[bed_name] + list(bed_name_other)].iloc[bed_review_idx[bed_name]], p], axis=1).sort_values(by=bed_name, ascending=False)

keys = list(beds.keys())