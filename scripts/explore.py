import pandas as pd
import numpy as np
import gzip
import os


# Carga de los archivos
PROJECT_ROOT = os.getcwd()
DATA_PATH = os.path.join(PROJECT_ROOT, "data", "raw")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "processed")

# Protein INFO
protein_info_csv = pd.read_csv(os.path.join(DATA_PATH, "9606.protein.info.v12.5.txt.gz"), sep="\t")
protein_info_csv.head()
protein_info_csv.columns
protein_info_csv.sample(5)
protein_info_csv.dtypes
protein_info_csv.shape
protein_info_csv.info()

# Tenemos 19698 proteinas, sin valores nulos
# Tenemos 4 columnas informativas: 1) STR id de proteina, parece que esta constituida por dos STR a su vez, 2) STR nombre de la proteina, 3) INT tamaño de la proteina, 4) STR informacion adicional

# Protein LINKS
protein_links_csv = pd.read_csv(os.path.join(DATA_PATH, "9606.protein.links.v12.5.txt.gz"), sep=r"\s+")
protein_links_csv.head()
protein_links_csv.columns
protein_links_csv.sample(5)
protein_links_csv.dtypes
protein_links_csv.shape
protein_links_csv.info()

# Tenemos 10826054 interacciones descritas
# Tenemos 3 columnas: 1) y 2) STR proteinas de la interaccion, 3) INT score de la interaccion

# Protein aliases
protein_alias_csv = pd.read_csv(os.path.join(DATA_PATH, "9606.protein.aliases.v12.5.txt.gz"), sep="\t")
protein_alias_csv.head()
protein_alias_csv.columns
protein_alias_csv.sample(5)
protein_alias_csv.dtypes
protein_alias_csv.shape
protein_alias_csv.info()
protein_alias_csv[protein_alias_csv["source"].str.contains("UniProt")]["source"].unique()
protein_alias_csv[protein_alias_csv["source"] == "UniProt_AC"]

# Este csv nos va a ayudar a conectar los archivos de informacion sobre las proteinas con el csv de los farmacos.
# El csv de los farmacos, como veremos mas adelante, utiliza la codificacion de Uniprot AC, y los otros dos archivos de proteinas usan otra codificacion


# Uniprot mapping
uniprot_map = pd.read_csv(
	os.path.join(DATA_PATH, "chembl_uniprot_mapping.txt"),
	sep="\t",
	comment="#",
	header=None,
	names=["uniprot_accession", "chembl_id", "target_name", "target_type"],
)
uniprot_map.head()
uniprot_map.columns
uniprot_map.dtypes
uniprot_map.sample(5)

# Contiene informacion sobre que proetine se ve target de que farmaco.

# Transformaciones de los CSV
# Del csv de alias, nos tenemos que quedar solo con los alias que nos interesan, que es UNIPROT AC
# Tenemos que filtrar del csv de links de manera que nos quedamos solo con los scores mas altos (>=700)
# Merge de protein links con los alias
# Merge con protein info csv
# uniprot_1  uniprot_2  combined_score  gene_name_1  gene_name_2
# P12345     Q67890     850             TP53         BRCA1
# P11111     Q22222     920             APP          PSEN1

# Filtro de UNIPROT AC
protein_alias_csv_filtered = protein_alias_csv[protein_alias_csv["source"] == "UniProt_AC"].drop_duplicates(subset="#string_protein_id", keep="first")
protein_alias_csv_filtered["source"].unique()
protein_alias_csv_filtered.head()

# Filtro combined_score
protein_links_csv_filtered = protein_links_csv[protein_links_csv["combined_score"] >= 700]
protein_links_csv_filtered.sample(10)

# Combinamos csvs de alias y links
protein_links_alias = protein_links_csv_filtered.merge(protein_alias_csv_filtered, left_on="protein1", right_on="#string_protein_id")
protein_links_alias = protein_links_alias.merge(protein_alias_csv_filtered, left_on="protein2", right_on="#string_protein_id")
protein_links_alias = protein_links_alias.rename(columns={"alias_x": "uniprot_1", "alias_y": "uniprot_2"})
protein_links_alias = protein_links_alias.drop(columns=["#string_protein_id_x", "source_x", "source_y", "#string_protein_id_y"])

# merge con info para saber el nombre
complete_protein_info = protein_links_alias.merge(
    protein_info_csv[["#string_protein_id", "preferred_name"]],
    left_on="protein1",
    right_on="#string_protein_id"
).merge(
    protein_info_csv[["#string_protein_id", "preferred_name"]],
    left_on="protein2",
    right_on="#string_protein_id",
    suffixes=("_1", "_2")
)
complete_protein_info = complete_protein_info.drop(columns=["#string_protein_id_1", "#string_protein_id_2", "protein1", "protein2"])