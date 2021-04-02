import pandas as pd
import synapseclient
import json


def get_data(syn, synapse_id):
    synapse_entity = syn.get(synapse_id)
    while isinstance(synapse_entity, synapseclient.entity.Link):
        synapse_entity = syn.get(synapse_entity.linksTo['targetId'])

    return pd.read_csv(synapse_entity.path, sep='\t', comment='#', low_memory=False)


def genie_fusion_counts(synapse_credentials_file, synapse_genie_release_id, filename=None):
    # synapse client
    syn = synapseclient.Synapse()
    # login, credential file: {"email": "", "password": ""}
    syn.login(**json.loads(open(synapse_credentials_file, 'r').read()))

    # make dict of entities in genie release, keyed by 'name'
    # provide syn id of folder of the release to be examined
    genie_syn_entities = {entity['name']: entity for entity in syn.getChildren(synapse_genie_release_id)}
    # genie_syn_entities.keys()

    # get fusion data
    fusions_df = get_data(syn, genie_syn_entities['data_fusions.txt'])
    # get sample data
    sample_df = get_data(syn, genie_syn_entities['data_clinical_sample.txt'])
    # seq assay id info
    assay_df = get_data(syn, genie_syn_entities['assay_information.txt'])

    # all types of alterations
    # alt_types = np.unique([_ for l in seq_df['alteration_types'].str.split(';').values for _ in l])
    alt_type = 'structural_variants'

    # list of seq assay types with reported fusion data
    fusion_samples = pd.merge(fusions_df['Tumor_Sample_Barcode'].drop_duplicates(),
                              sample_df[['SAMPLE_ID', 'SEQ_ASSAY_ID']],
                              left_on='Tumor_Sample_Barcode',
                              right_on='SAMPLE_ID',
                              how='left')
    fusion_samples = fusion_samples['SEQ_ASSAY_ID'].value_counts().reset_index()
    fusion_samples.rename(columns={'index': 'SEQ_ASSAY_ID', 'SEQ_ASSAY_ID': 'Samples with fusion calls'}, inplace=True)
    # add flag as to whether SEQ_ASSAY_ID is described as covering the alt_type
    fusion_samples['Alteration type includes "structural_variants"'] = fusion_samples['SEQ_ASSAY_ID'].isin(assay_df['SEQ_ASSAY_ID'][assay_df['alteration_types'].str.contains(alt_type)])

    # write to file
    if filename is not None:
        fusion_samples.to_excel(filename, index=False)


if __name__ == "__main__":
    # execute only if run as a script
    genie_fusion_counts('.synapse_credentials', 'syn23775220', 'fusion_panels_9.6.xlsx')
