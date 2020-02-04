import numpy as np
import pandas as pd
from tqdm import tqdm
import pickle
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import pyranges as pr


# load maf
# columns to load from maf
# cols = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 't_alt_count', 't_ref_count']
# maf = pd.read_csv('data/data_mutations_extended_7.6-consortium.txt', sep='\t', usecols=cols, low_memory=False)
maf = pd.read_csv('data/data_mutations_extended_7.6-consortium.txt', sep='\t', low_memory=False)
# idx = maf.duplicated(keep=False)
# t = maf.loc[idx]

# load sample table
sample = pd.read_csv('data/data_clinical_sample_7.6-consortium.txt', sep='\t', skiprows=4)

# load bed file
gi = pd.read_csv('data/genomic_information_7.6-consortium.txt', sep='\t', low_memory=False)
beds = dict(list(gi[['Chromosome', 'Start_Position', 'End_Position', 'SEQ_ASSAY_ID']].groupby('SEQ_ASSAY_ID')))

# add assay id to amf
maf = pd.merge(maf, sample[['SAMPLE_ID', 'SEQ_ASSAY_ID']], how='left', left_on='Tumor_Sample_Barcode', right_on='SAMPLE_ID')
# make unique variant string
maf['var_str'] = maf['Chromosome'].astype(str) + '_' + maf['Start_Position'].astype(str) + '_' + maf['End_Position'].astype(str) + '_' + maf['Reference_Allele'].astype(str) + '_' + maf['Tumor_Seq_Allele2'].astype(str)

# TODO: Need to update to reflect GENIE germline filtering approach
idx_germline = (maf['FILTER'].str.contains('common_variant').fillna(False)) | (maf.iloc[:, -10:-1].max(axis=1, skipna=True) > 0.001)

# counts of mutations across assay ID (nan = not in data)
variant_counts = maf[~idx_germline].groupby(['var_str', 'SEQ_ASSAY_ID', 'Tumor_Sample_Barcode']).size().to_frame().reset_index()
sample_counts = variant_counts.groupby(['var_str', 'SEQ_ASSAY_ID']).size().to_frame().reset_index().pivot(index='var_str', columns='SEQ_ASSAY_ID', values=0)

# get counts of sample per assay
panel_counts = sample.groupby('SEQ_ASSAY_ID').size()

# get unique mutations
var_uniq = maf.loc[~idx_germline, ['var_str', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']].drop_duplicates().set_index('var_str')
# make PyRange for unique variants, include index (var_str) for setting value in orginal df
var_pr = pr.PyRanges(var_uniq.reset_index()[['Chromosome', 'Start_Position', 'End_Position', 'var_str']].rename(columns={'Start_Position': 'Start', 'End_Position': 'End'}))

# find overlap of variants in various bed files and add to unique variant table
for bed_name in beds.keys():
    var_uniq[bed_name] = False
    var_uniq.loc[var_pr.overlap(pr.PyRanges(chromosomes=beds[bed_name]['Chromosome'], starts=beds[bed_name]['Start_Position'], ends=beds[bed_name]['End_Position'])).var_str, bed_name] = True

# is nan to true zeros based on coverage of panel for specific variant
sample_counts.values[(sample_counts.isna() & var_uniq.iloc[:, 6:]).values] = 0

bed_review_idx = dict()
# criteria for flagging variants for assay bias
for bed_name in beds.keys():
    bed_name_other = np.setdiff1d(list(beds.keys()), bed_name)
    bed_review_idx[bed_name] = np.where((sample_counts[bed_name] >= 10) & (sample_counts[bed_name_other].sum(axis=1, skipna=True) < 10) & ((~sample_counts[bed_name_other].isna()).sum(axis=1) >= 3))[0]

bed_name = 'MSK-IMPACT468'
bed_name = 'JHU-50GP'
bed_name = 'UHN-555-V1'
bed_name = 'PHS-FOCUS-V1'

bed_name_other = np.setdiff1d(list(beds.keys()), bed_name)
p = sample_counts.iloc[bed_review_idx[bed_name]].apply(lambda x: pd.Series(fisher_exact([[x[bed_name], panel_counts[bed_name]], [x[bed_name_other].sum(skipna=True), panel_counts[bed_name_other[(~x[bed_name_other].isna())]].sum()]], alternative='greater')[1], index=['p-value']), axis=1)
review_table = pd.concat([var_uniq.loc[sample_counts.index[bed_review_idx[bed_name]]].iloc[:, :6], sample_counts.iloc[bed_review_idx[bed_name]]], axis=1).sort_values(by=bed_name, ascending=False)
review_table['n_' + bed_name] = panel_counts[bed_name]
review_table['not_' + bed_name] = review_table[bed_name_other].sum(axis=1, skipna=True)
review_table['n_not_' + bed_name] = (~review_table[bed_name_other].isna()).apply(lambda x: panel_counts[x.index[x].values].sum(), axis=1)
review_table['freq_' + bed_name] = review_table[bed_name] / review_table['n_' + bed_name]
review_table['freq_not' + bed_name] = review_table['not_' + bed_name] / review_table['n_not_' + bed_name]
review_table['fisher_p'] = review_table.apply(lambda x: fisher_exact([[x[bed_name], x['n_' + bed_name] - x[bed_name]], [x['not_' + bed_name], x['n_not_' + bed_name] - x['not_' + bed_name]]], alternative='greater')[1], axis=1)
review_table = review_table[list(review_table.columns[:6]) + [bed_name] + list(review_table.columns[-6:]) + list(bed_name_other)]

review_table.to_csv(bed_name + '_for_review.tsv', sep='\t', index_label='index')

