import pandas
from pathlib import Path
from tools.mapping.ensembl_uniprot import idmapping_uniprot_mapping

egfr_directory = Path(__file__).parent.resolve() / '..'

def main():
    # Re-read the uniprot nodes from `process_gold_standard.py`
    nodes = (egfr_directory / 'processed' / 'gold-standard-nodes-uniprot.txt').read_text().splitlines()
    # and the prizes from `process_prizes.py`
    prizes = pandas.read_csv(egfr_directory / 'processed' / 'prizes-uniprot.txt', sep='\t')

    # We grab our UniProt <-> ENSP mapping
    idmapping_df = idmapping_uniprot_mapping(egfr_directory / 'raw' / 'HUMAN_9606_idmapping_selected.tsv')

    # and map the nodes
    idmapping_nodes_df = pandas.DataFrame(nodes, columns=['UniProtKB-ID']).merge(idmapping_df, on='UniProtKB-ID', how='left')
    idmapping_nodes_df = idmapping_nodes_df.drop(columns=['UniProtKB-ID', 'UniProtKB-AC', 'Ensembl'])
    idmapping_nodes_df = idmapping_nodes_df[~idmapping_nodes_df['Ensembl_PRO'].isna()]
    nodes = idmapping_nodes_df['Ensembl_PRO'].astype(str).to_list()
    (egfr_directory / 'processed' / 'gold-standard-nodes.txt').write_text("\n".join(nodes))

    # and the prizes
    idmapping_prizes_df = prizes.merge(idmapping_df, left_on='NODEID', right_on="UniProtKB-ID", how='inner')
    idmapping_prizes_df = idmapping_prizes_df.drop(columns=['UniProtKB-ID', 'UniProtKB-AC', 'Ensembl', 'NODEID'])
    idmapping_prizes_df = idmapping_prizes_df[~idmapping_prizes_df['Ensembl_PRO'].isna()]
    idmapping_prizes_df = idmapping_prizes_df.rename(columns={'Ensembl_PRO': 'NODEID'})
    idmapping_prizes_df = idmapping_prizes_df[["NODEID", "prize", "active", "dummy", "source"]]
    idmapping_prizes_df.to_csv(egfr_directory / 'processed' / 'prizes.txt', sep='\t', index=False)

if __name__ == "__main__":
    main()