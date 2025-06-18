import pickle

with open('Pickles/UniprotIDs.pkl', 'rb') as file:
        UniprotIDs = pickle.load(file)

UIDs = UniprotIDs["UniprotIDs"]
UMap= UniprotIDs["UniprotMap"]

with open('Pickles/NodeIDs.pkl','rb') as file2:
        prizes = pickle.load(file2)

prize_05 = prizes["prize_05"]
prize_060 = prizes["prize_060"]

prize_05['Uniprot'] = prize_05['Uniprot'].apply(lambda x: UMap.get(x))
prize_060['Uniprot'] = prize_060['Uniprot'].apply(lambda x: UMap.get(x))

prize_05.columns = ['NODEID','prize']
prize_060.columns = ['NODEID','prize']

df = {
        "prize_05":prize_05,
        "prize_060":prize_060
}

with open("Pickles/SPRAS_Data.pkl","wb") as file:
        pickle.dump(df,file)