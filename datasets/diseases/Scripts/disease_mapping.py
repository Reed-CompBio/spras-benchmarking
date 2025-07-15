import pandas as pd 
import requests
import pickle
import sys

def main():

    # Set pandas df display options
    # pd.set_option('display.max_columns',None)  
    pd.set_option('display.max_rows', None)   

    # Read in raw data files
    text_mining = pd.read_csv('datasets/diseases/raw/human_disease_textmining_filtered.tsv',sep = '\t')
    text_mining.columns= ['geneID','geneName','diseaseID','diseaseName','zScore','confidenceScore','sourceUrl']

    knowlege = pd.read_csv('datasets/diseases/raw/human_disease_knowledge_filtered.tsv',sep = '\t')
    knowlege.columns= ['geneID','geneName','diseaseID','diseaseName','sourceDB','evidenceType','confidenceScore']

    experiments = pd.read_csv('datasets/diseases/raw/human_disease_experiments_filtered.tsv',sep = '\t')
    experiments.columns= ['geneID','geneName','diseaseID','diseaseName','sourceDB','sourceScore','confidenceScore']

    integrated =  pd.read_csv('datasets/diseases/raw/human_disease_integrated_full.tsv',sep = '\t')
    integrated.columns = ['geneID','geneName','diseaseID','diseaseName','integrated_confidenceScore']
    integrated = integrated.drop(columns=['geneName','diseaseName'])

    # Crop dataframes to remove unused columns and remove duplicate values
    tiga = pd.read_csv('datasets/diseases/raw/tiga_gene-trait_stats.tsv',sep='\t')
    tiga_cropped = tiga[['ensemblId','efoId','trait','n_snp','n_snpw']]
    # duplicate_mask = tiga_cropped.duplicated(subset=['ensemblId','trait'])
    # tiga_duplicates = tiga_cropped[duplicate_mask]
    # print(tiga_duplicates[['ensemblId','trait']])
    tiga_cropped = tiga_cropped.drop_duplicates(subset=['ensemblId','trait'])

    human_do = pd.read_csv('datasets/diseases/raw/HumanDO.tsv',sep = '\t')
    human_do = human_do.drop_duplicates(subset='label')
    do_cropped = human_do[['id','label']]

    # Merge tiga dataframe with human disease ontology data to assign each tiga trait to a disease ontology code
    tiga_do = tiga_cropped.merge(do_cropped,left_on='trait',right_on='label',how='inner',validate='m:1')

    #Map geneIDs in experiments from ENSP to ENSG format
    r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism':'hsapiens',
        'target':'ENSG',
        'query': list(experiments['geneID']),
        }
    )
    results = r.json()['result']
    mapping ={}
    for x in results:
        if x['converted'] != 'None':
            mapping.update({x['incoming']:x['converted']})

    mapping_df = pd.DataFrame.from_dict(mapping.items())
    mapping_df.columns = ['ENSP','ENSG']
    experiments_mapped = experiments.merge(mapping_df,left_on='geneID',right_on='ENSP',how ='inner',validate='m:1')

    # Merge tiga data with experimental data using both gene id and disease id as keys
    tiga_experiments = experiments_mapped.merge(tiga_do,left_on=['ENSG','diseaseID'],right_on=['ensemblId','id'],how='inner',validate='1:1')
    tiga_experiments = tiga_experiments[['ENSP','ENSG','geneName','trait','efoId','diseaseID','sourceScore','confidenceScore','n_snp','n_snpw']]
    # print(tiga_experiments[tiga_experiments.isna().any(axis=1)])
    # sys.exit()

    # Get string ID's for each protein in tiga_experiments
    string_api_url = "https://version-12-0.string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"
    str_params = {
        "identifiers" : "\r".join(list(tiga_experiments['ENSP'])), 
        "species" : 9606, 
        "echo_query" : 1, 
    }
    request_url = "/".join([string_api_url, output_format, method])
    string_results = requests.post(request_url, data=str_params)

    string_map = {}
    for line in string_results.text.strip().split("\n"):
        l = line.split("\t")
        string_map.update({l[0]:l[2]})
    string_df = pd.DataFrame.from_dict(string_map.items())
    string_df.columns = ['ENSP','str_id']

    # Merge string ID's with tiga_experiments such that each protein has a corresponding string ID
    merged_df = tiga_experiments.merge(string_df,on = 'ENSP',how ='inner')
    final_df = merged_df.merge(integrated,left_on=['ENSP','diseaseID'],right_on=['geneID','diseaseID'],how = 'inner',validate='1:1')
    # final_df = final_df.drop(columns='geneName_y')

    # print(final_df)

    #Convert final df to dictionary where each key is a trait 
    trait_group = final_df.groupby('trait')
    trait_dict = {k:v for k,v in trait_group}

    df = {
        "final_df": final_df,
        "trait_dict": trait_dict
    }

    with open("datasets/diseases/viz/pklDat.pkl","wb") as file:
        pickle.dump(df,file)
   
    return 


if __name__ == '__main__':
    main()