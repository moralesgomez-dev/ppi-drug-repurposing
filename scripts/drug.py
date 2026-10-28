import pandas as pd
import os
from pyvis.network import Network

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

# Extraemos los vecinos de primer grado de los genes Alzheimer que no son genes Alzheimer.
alzheimer_ids = set(genes_in_network["uniprot_id"].dropna())
print(f"Proteínas Alzheimer únicas en la red: {len(alzheimer_ids)}")
neighbor_ids = {
    neighbor
    for alzheimer_id in alzheimer_ids
    for neighbor in G_giant.neighbors(alzheimer_id)
    if neighbor not in alzheimer_ids
}
neighbor_genes = pd.DataFrame({"uniprot_id": sorted(neighbor_ids)})
neighbor_genes = neighbor_genes.merge(degree_df, on="uniprot_id", how="left")
print(f"Proteínas vecinas no Alzheimer: {len(neighbor_genes)}")

# Cruzamos los vecinos con el mapeo de ChEMBL a UniProt generado en explore.py.
uniprot_map = pd.read_csv(os.path.join(DATA_PATH, "processed\\uniprot_map.csv"))
uniprot_map.head()

# Cruzamos las proteínas vecinas con uniprot_map para añadir la columna de chembl_id.
neighbors_with_targets = neighbor_genes.merge(
    uniprot_map,
    left_on="uniprot_id",
    right_on="uniprot_accession",
    how="left",
)
neighbors_with_targets.head()

# Vemos presecia de NaN en chembl_id
print("Número de filas con NaN en chembl_id:", neighbors_with_targets["chembl_id"].isna().sum())
# Estas proteínas vecinas no tienen una diana asociada en la base de datos de ChEMBL.

# Nos quedamos con las proteínas que tienen una diana asociada en ChEMBL.
repurposing_candidates = neighbors_with_targets.dropna(subset=["chembl_id"])

# Conservamos solo las dianas SINGLE PROTEIN y ordenamos por centralidad.
repurposing_candidates = repurposing_candidates[
    repurposing_candidates["target_type"] == "SINGLE PROTEIN"
]
repurposing_candidates = repurposing_candidates.drop(columns=["uniprot_accession"])
repurposing_candidates = repurposing_candidates.sort_values(
    by="degree_centrality",
    ascending=False,
)
print(f"Registros de dianas ChEMBL candidatos: {len(repurposing_candidates)}")
print(
    "Proteínas vecinas únicas con diana ChEMBL: "
    f"{repurposing_candidates['uniprot_id'].nunique()}"
)

# Guardamos el resultado final en un csv
repurposing_candidates.to_csv(
    os.path.join(DATA_PATH, "processed\\alzheimer_neighbor_targets_with_chembl.csv"),
    index=False,
)

# Visualizamos las 30 dianas candidatas con mayor centralidad y sus conexiones
# directas con proteínas asociadas al Alzheimer.
TOP_TARGET_COUNT = 30
top_targets = repurposing_candidates.drop_duplicates(
    subset=["uniprot_id"]
).head(TOP_TARGET_COUNT)
top_target_ids = set(top_targets["uniprot_id"])
candidate_edges = [
    (source, target, edge_data)
    for source, target, edge_data in G_giant.edges(data=True)
    if (
        source in alzheimer_ids
        and target in top_target_ids
    )
    or (
        target in alzheimer_ids
        and source in top_target_ids
    )
]
connected_alzheimer_ids = {
    node
    for source, target, _ in candidate_edges
    for node in (source, target)
    if node in alzheimer_ids
}

network = Network(
    height="800px",
    width="100%",
    directed=False,
    notebook=False,
    bgcolor="#ffffff",
    font_color="#222222",
)

gene_symbols = (
    genes_in_network.dropna(subset=["uniprot_id"])
    .groupby("uniprot_id")["symbol"]
    .agg(lambda symbols: ", ".join(sorted(set(symbols))))
    .to_dict()
)
for alzheimer_id in connected_alzheimer_ids:
    symbol = gene_symbols.get(alzheimer_id, alzheimer_id)
    network.add_node(
        alzheimer_id,
        label=symbol,
        title=f"Proteína asociada al Alzheimer | UniProt: {alzheimer_id}",
        color="#3182bd",
        shape="dot",
    )

for _, row in top_targets.iterrows():
    protein_id = row["uniprot_id"]
    target_details = (
        f"Proteína vecina con diana ChEMBL | "
        f"UniProt: {protein_id} | "
        f"Diana: {row['target_name']} | "
        f"ChEMBL: {row['chembl_id']} | "
        f"Centralidad de grado: {row['degree_centrality']:.6f}"
    )
    network.add_node(
        protein_id,
        label=protein_id,
        title=target_details,
        color="#2ca25f",
        shape="diamond",
    )

for source, target, edge_data in candidate_edges:
    network.add_edge(
        source,
        target,
        title=f"Puntuación combinada: {edge_data.get('combined_score', 'N/D')}",
    )

network.show_buttons(filter_=["physics"])
FIGURES_PATH = os.path.join(PROJECT_ROOT, "results", "figures")
os.makedirs(FIGURES_PATH, exist_ok=True)
network_path = os.path.join(
    FIGURES_PATH,
    "alzheimer_neighbor_targets.html",
)
network.write_html(network_path, open_browser=False, notebook=False)
print(f"Red de genes Alzheimer y dianas vecinas guardada en: {network_path}")
