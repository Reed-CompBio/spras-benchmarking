import pickle 
import pandas as pd

def main():
    

    with open('datasets/diseases/Pickles/GS.pkl', 'rb') as file:
        pklGS = pickle.load(file)

    with open('datasets/diseases/Pickles/Inputs.pkl', 'rb') as file:
        pklInput = pickle.load(file)

    
    GS_string_df = pklGS['GS_string_df']
    tiga_string_df = pklInput['tiga_string_df']


    GS_string_df = GS_string_df[GS_string_df['diseaseID'].isin(tiga_string_df['id'])]
    GS_combined_group = GS_string_df.groupby('diseaseName')
    GS_combined_dict = {k:v for k,v in GS_combined_group}


    tiga_filtered = tiga_string_df[tiga_string_df['id'].isin(GS_string_df['diseaseID'])]
    tiga_group = tiga_filtered.groupby('trait')
    tiga_dict = {k:v for k,v in tiga_group}
    tiga_count = {x:len(tiga_dict[x]) for x in tiga_dict.keys()}
    tiga_count_threshold = {k:v  for (k,v) in tiga_count.items() if (v > 10)}

    tiga_threshold = tiga_filtered.loc[tiga_filtered['trait'].isin(list(tiga_count_threshold.keys()))]

    tiga_prizes = tiga_threshold.groupby('trait')
    tiga_prize_dict = {k:v for k,v in tiga_prizes}

    for disease in tiga_prize_dict.keys():
        df = tiga_prize_dict[disease]
        df = df[['str_id','n_snpw']]
        df = df.rename(columns={'str_id':'NODEID','n_snpw':'prize'})
        df.to_csv(f"datasets/diseases/prize_files/{disease.replace(' ','_')}_prizes.txt",sep = '	',index=False)

    for disease in GS_combined_dict.keys():
        df = GS_combined_dict[disease]
        df = df[['str_id']]
        df.to_csv(f"datasets/diseases/GS_files/{disease.replace(' ','_')}_GS.txt",sep = '	',index=False,header=None)

    string = pd.read_csv('datasets/diseases/raw/9606.protein.links.v12.0.txt',sep = ' ',skiprows=[0],header=None)
    string = string[string.iloc[:,2]>900]
    string = string.iloc[:,[0,1]]
    string[len(string.columns)] = 1
    string.to_csv('datasets/diseases/raw/string_interactome.txt',sep = '\t',index=False,header=None)

    return



if __name__ == '__main__':
    main()