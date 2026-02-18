import pandas
from pathlib import Path
from tools.mapping.ensembl_uniprot import idmapping_uniprot_mapping

egfr_directory = Path(__file__).parent.resolve() / '..'

def main():
    # First, we remove all PSUEDONODES (and any duplicates)
    nodes = (egfr_directory / 'raw' / 'eight-egfr-reference-all.txt').read_text().splitlines()
    nodes = list(set([node for node in nodes if not node.endswith("_PSEUDONODE")]))

    # Then, we map our UniProt nodes to ENSP.
    idmapping_df = idmapping_uniprot_mapping(egfr_directory / 'raw' / 'HUMAN_9606_idmapping_selected.tsv')
    idmapping_df = pandas.DataFrame(nodes, columns=['UniProtKB-ID']).merge(idmapping_df, on='UniProtKB-ID', how='left')
    idmapping_df = idmapping_df.drop(columns=['UniProtKB-ID', 'UniProtKB-AC', 'Ensembl'])
    idmapping_df = idmapping_df[~idmapping_df['Ensembl_PRO'].isna()]
    nodes = idmapping_df['Ensembl_PRO'].astype(str).to_list()

    (egfr_directory / 'processed').mkdir(exist_ok=True)
    (egfr_directory / 'processed' / 'gold-standard-nodes.txt').write_text("\n".join(nodes))

if __name__ == "__main__":
    main()
