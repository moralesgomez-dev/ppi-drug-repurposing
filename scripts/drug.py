import pandas as pd
import os

from scripts.grafo import G_giant

PROJECT_ROOT = os.getcwd()
DATA_PATH = os.path.join(PROJECT_ROOT, "data")  

# Carga del archivo
genes = pd.read_csv(os.path.join(DATA_PATH, "raw\\OT-MONDO_0004975-associated-targets-26_10_2026-v26_09 (2).tsv"), sep='\t')

# Vemos la estructura que tiene el archivo descargado de OpenTargets
print(genes.head())
print(genes.columns)

# Solo nos interesan las columnas de symbol, globalScore, GWASCredibleSets y uniprotLiterature
genes = genes[['symbol', 'globalScore', 'gwasCredibleSets', 'uniprotLiterature']]
print(genes.head())

# En este caso tenemos el nombre de los genes (symbol), pero no tenemos el identificador de Uniprot, que es el que necesitamos para poder cruzar con la información de proteínas. 
# Por lo tanto, vamos a hacer un mapeo de los símbolos de genes a sus identificadores de Uniprot.
# Tenemos ese mapeo en protein_info_csv, que es el archivo que hemos generado en el script grafo.py. 
# Por lo tanto, vamos a cargar ese archivo y hacer un merge con el dataframe genes para obtener los identificadores de Uniprot.
protein_info = pd.read_csv(os.path.join(DATA_PATH, "processed\\protein_info_csv"))
protein_info.head()



# Crear una tabla gen -> UniProt usando ambos extremos de cada interacción
mapping_1 = protein_info[["gene_1", "uniprot_1"]].rename(
    columns={"gene_1": "symbol", "uniprot_1": "uniprot_id"}
)
mapping_2 = protein_info[["gene_2", "uniprot_2"]].rename(
    columns={"gene_2": "symbol", "uniprot_2": "uniprot_id"}
)

gene_uniprot = (
    pd.concat([mapping_1, mapping_2], ignore_index=True)
    .dropna(subset=["symbol", "uniprot_id"])
    .drop_duplicates(subset=["symbol", "uniprot_id"])
)

# Revisar si algún símbolo tiene más de un UniProt ID
mapping_counts = gene_uniprot.groupby("symbol")["uniprot_id"].nunique()
ambiguous = mapping_counts[mapping_counts > 1]
print("Símbolos con varios UniProt IDs:")
print(ambiguous)

# Añadir el UniProt ID a los genes de Open Targets
genes_mapped = genes.merge(gene_uniprot, on="symbol", how="left")
genes_mapped.head()

#TODO: limpiar genes_mapped, cruzar con nodos del grado y añadir degree, añadir chembl farmcaos y limpiar todo y tener una tabla final

genes_mapped = genes_mapped[["symbol", "uniprot_id", "globalScore"]]
genes_mapped.head()

# Cargamos el df de degree centrality generado en grafo.py
degree_df = pd.read_csv(os.path.join(DATA_PATH, "processed\\degree_centrality.csv"))
degree_df.head()

# Lo cruzamos con genes_mapped para añadir la columna de degree_centrality
genes_mapped = genes_mapped.merge(degree_df, on="uniprot_id", how="left")
genes_mapped.head()
genes_mapped.sample(5)

# Vemos presecia de NaN en degree_centrality
print("Número de filas con NaN en degree_centrality:", genes_mapped["degree_centrality"].isna().sum())
# Esto nos indica que hay genes que no están presentes en el grafo de proteínas, por lo tanto no tienen un valor de centralidad de grado.
genes_in_network = genes_mapped.dropna(subset=["degree_centrality"])
print(f"Genes Alzheimer en la red PPI: {len(genes_in_network)}")

# Por ultimo falta cruzarlo con chEMBL para ver que farmacos estan asociados a estos genes. 
# Cargamos el archivo de mapeo de ChEMBL a UniProt que generamos en explore.py
uniprot_map = pd.read_csv(os.path.join(DATA_PATH, "processed\\uniprot_map.csv"))
uniprot_map.head()

# Cruzamos genes_in_network con uniprot_map para añadir la columna de chembl_id
genes_with_drugs = genes_in_network.merge(uniprot_map, left_on="uniprot_id", right_on="uniprot_accession", how="left")
genes_with_drugs.head()
genes_with_drugs.sample(5)

# Vemos presecia de NaN en chembl_id
print("Número de filas con NaN en chembl_id:", genes_with_drugs["chembl_id"].isna().sum())
# Esto nos indica que hay genes que no tienen un fármaco asociado en la base de datos de ChEMBL. Por lo tanto, no todos los genes de Alzheimer tienen un fármaco asociado.

# Nos quedamos solo con los genes que tienen un fármaco asociado y los ordanamos por globalScore y degree_centrality
final_genes = genes_with_drugs.dropna(subset=["chembl_id"])
final_genes = final_genes.sort_values(by=["globalScore"], ascending=[False])
final_genes.head()

# Ademas nos quedamos con SINGLE PROTEIN targets, que son los que nos interesan para el repurposing de fármacos. Quitamos uniprot_acession
final_genes = final_genes[final_genes["target_type"] == "SINGLE PROTEIN"]
final_genes = final_genes.drop(columns=["uniprot_accession"])
final_genes.head()

# Guardamos el resultado final en un csv
final_genes.to_csv(os.path.join(DATA_PATH, "processed\\alzheimer_genes_with_drugs.csv"), index=False)



