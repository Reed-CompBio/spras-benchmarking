import pandas as pd
from pathlib import Path
import sys

dir_path = Path(__file__).parent.resolve()
diseases_path = Path(dir_path, "..")

CONFIDENCE_SCORE_MINIMUM = 4
GENE_SET_SIZE_MINIMUM = 10

'''
Generates gold standard disease-gene associations from the DISEASES
database using the text mining and knowledge channels. It relies on 
a confidence score threshold to establish high confidence disease gene
pairs and then a size threshold to establish the diseases with enough
high confidence disease-gene pairs.
'''
def main():
    print('Reading data...')
    # Get our data from `fetch.py`
    text_mining = pd.read_csv(diseases_path / "raw" / "human_disease_textmining_full.tsv", sep="\t")
    knowledge = pd.read_csv(diseases_path / "raw" / "human_disease_knowledge_full.tsv", sep="\t")
    # and correctly map their columns
    text_mining.columns = ["geneID", "geneName", "diseaseID", "diseaseName", "zScore", "confidenceScore", "sourceUrl"]
    knowledge.columns = ["geneID", "geneName", "diseaseID", "diseaseName", "sourceDB", "evidenceType", "confidenceScore"]
    print(f' {len(text_mining)} text mining and {len(knowledge)} knowledge associations.')

    # retain the columns needed
    text_mining = text_mining[["geneID", "geneName", "diseaseID", "diseaseName", "confidenceScore"]]
    knowledge = knowledge[["geneID", "geneName", "diseaseID", "diseaseName", "confidenceScore"]]

    # drop disease-gene associations if their confidence scores 
    # are less than CONFIDENCE_SCORE_MINIMUM
    # doing this here so the dataframes are substantially smaller.
    text_mining = text_mining[text_mining["confidenceScore"] >= CONFIDENCE_SCORE_MINIMUM]
    knowledge = knowledge[knowledge["confidenceScore"] >= CONFIDENCE_SCORE_MINIMUM]
    print(f' {len(text_mining)} text mining and {len(knowledge)} knowledge associations with confidenceScore >= {CONFIDENCE_SCORE_MINIMUM}.')

    print('\nCombining Text Mining and Knowledge Channels...')
    # Take the maximum of the text mining and knowledge confidence score as the final confidence score.
    inner = text_mining.merge(knowledge, on=["geneID", "diseaseID"], how="inner")
    inner["confidenceScore"] = inner.apply(lambda x: max(x.confidenceScore_x, x.confidenceScore_y), axis=1)
    inner = inner.rename(columns={"geneID_x": "geneID", "geneName_x": "geneName", "diseaseName_x": "diseaseName"})
    inner = inner[["geneID","geneName","diseaseID","diseaseName","confidenceScore"]]
    print(f' {len(inner)} associations found in both text mining and knowledge channels. Maximum score retained.')

    # Take the text mining score for any disease-gene association that does not have a knowledge score.
    txt_only = text_mining[~text_mining["geneID"].isin(knowledge["geneID"])]
    print(f' {len(txt_only)} associations from text mining only.')

    # Take the knowledge score for any disease-gene association that does not have a text mining score.
    kn_only = knowledge[~knowledge["geneID"].isin(text_mining["geneID"])]
    print(f' {len(kn_only)} associations from knowledge only.')

    # combine all genes that have at least one text mining or knowledge score.
    GS_high_conf = pd.concat([inner, txt_only, kn_only])
    print(f' {len(GS_high_conf)} total high confidence disease-gene associations.')
    GS_high_conf.to_csv("test.tsv", index=False)

    print('\nFiltering diseases...')

    # Require that disease have at least GENE_SET_SIZE_MINIMUM high-confidence gene associations.
    GS_diseases = GS_high_conf.groupby("diseaseID")
    print(f' The high-confidence disease-gene pairs correspond to {len(GS_diseases)} distinct diseases.')
    GS_retained_diseases = set([k for k, v in GS_diseases if (len(v) >= GENE_SET_SIZE_MINIMUM)])
    print(f' There are {len(GS_retained_diseases)} diseases with at least {GENE_SET_SIZE_MINIMUM} high-confidence disease-gene pairs')
    
    # Filter the high-confidence disease-gene pairs to retain only the diseases that have been retained.
    GS_high_conf_retained_diseases = GS_high_conf[GS_high_conf["diseaseID"].isin(GS_retained_diseases)]
    print(f' There are {len(GS_high_conf_retained_diseases)} high-confidence disease-gene pairs from the {len(GS_retained_diseases)} diseases')

    ## TODO - do we need to filter by STRING??

    # Write output file gold_standard.csv
    GS_high_conf_retained_diseases.to_csv(diseases_path / "processed" / "gold_standard.csv", index=False)

if __name__ == "__main__":
    main()
