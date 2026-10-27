from html import escape
import os

import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from pyvis.network import Network

# Rutas del proyecto
PROJECT_ROOT = os.getcwd()
DATA_PATH = os.path.join(PROJECT_ROOT, "data")
PROCESSED_PATH = os.path.join(DATA_PATH, "processed")
RAW_PATH = os.path.join(DATA_PATH, "raw")
OUTPUT_PATH = os.path.join(
    PROJECT_ROOT, "results", "figures", "alzheimer_subgraph.html"
)
FIGURE_PATH = os.path.join(
    PROJECT_ROOT, "results", "figures", "alzheimer_top20_degree_centrality.png"
)

# Colores por tipo de nodo
COLORS = {
    "chembl_target": "#2ca25f",
    "no_chembl_target": "#3182bd",
    "hub": "#de2d26",
}

# Carga de los datos
genes = pd.read_csv(
    os.path.join(
        RAW_PATH, "OT-MONDO_0004975-associated-targets-26_10_2026-v26_09 (2).tsv"
    ),
    sep="\t",
    usecols=["symbol", "globalScore"],
)
protein_info = pd.read_csv(
    os.path.join(PROCESSED_PATH, "protein_info_csv"),
    usecols=[
        "uniprot_1",
        "uniprot_2",
        "gene_1",
        "gene_2",
        "combined_score",
    ],
)

degree_df = pd.read_csv(os.path.join(PROCESSED_PATH, "degree_centrality.csv"))
degree_df["uniprot_id"] = degree_df["uniprot_id"].astype(str)
chembl_targets = pd.read_csv(os.path.join(PROCESSED_PATH, "uniprot_map.csv"))

# Mapeo de genes Alzheimer a UniProt
giant_nodes = set(degree_df["uniprot_id"])
mapping_1 = protein_info[["gene_1", "uniprot_1"]].rename(
    columns={"gene_1": "symbol", "uniprot_1": "uniprot_id"}
)

mapping_2 = protein_info[["gene_2", "uniprot_2"]].rename(
    columns={"gene_2": "symbol", "uniprot_2": "uniprot_id"}
)

gene_mapping = (
    pd.concat([mapping_1, mapping_2], ignore_index=True)
    .dropna(subset=["symbol", "uniprot_id"])
    .drop_duplicates(subset=["symbol", "uniprot_id"])
)

alzheimer_nodes = (
    genes.merge(gene_mapping, on="symbol", how="inner")
    .drop_duplicates(subset=["symbol", "uniprot_id"])
)

alzheimer_nodes["uniprot_id"] = alzheimer_nodes["uniprot_id"].astype(str)
alzheimer_nodes = alzheimer_nodes[alzheimer_nodes["uniprot_id"].isin(giant_nodes)]
alzheimer_nodes = alzheimer_nodes.sort_values(["uniprot_id", "symbol"])

# Nos quedamos con el símbolo de gen más representativo para cada UniProt ID y el globalScore máximo
gene_summary = alzheimer_nodes.groupby("uniprot_id").agg(
    {"symbol": ", ".join, "globalScore": "max"}
)
gene_summary = gene_summary.reset_index()
if gene_summary["uniprot_id"].duplicated().any():
    raise ValueError("gene_summary contiene UniProt IDs duplicados")

# Asociación de dianas ChEMBL
chembl_targets = chembl_targets[
    chembl_targets["target_type"].eq("SINGLE PROTEIN")
    & chembl_targets["chembl_id"].notna()
]

chembl_targets = chembl_targets.drop_duplicates(
    subset=["uniprot_accession", "chembl_id", "target_name"]
)

chembl_targets["chembl_target"] = (
    chembl_targets["chembl_id"].astype(str)
    + " - "
    + chembl_targets["target_name"].astype(str)
)

chembl_by_uniprot = chembl_targets.groupby("uniprot_accession").agg(
    {"chembl_target": "; ".join}
)
chembl_by_uniprot = chembl_by_uniprot.reset_index()

# Merge de la información de dianas ChEMBL y degree centrality con el resumen de genes
gene_summary = gene_summary.merge(
    chembl_by_uniprot,
    left_on="uniprot_id",
    right_on="uniprot_accession",
    how="left",
)

gene_summary = gene_summary.merge(degree_df, on="uniprot_id", how="left")
gene_summary["chembl_target"] = gene_summary["chembl_target"].fillna("")

# Construcción del subgrafo y clasificación de nodos
top_50_hubs = degree_df.sort_values(
    "degree_centrality", ascending=False
).head(50)["uniprot_id"].tolist()
gene_summary["category"] = "no_chembl_target"

for index, row in gene_summary.iterrows():
    if row["uniprot_id"] in top_50_hubs:
        gene_summary.loc[index, "category"] = "hub"
    elif row["chembl_target"]:
        gene_summary.loc[index, "category"] = "chembl_target"

alzheimer_ids = gene_summary["uniprot_id"].tolist()
alzheimer_edges = protein_info[
    protein_info["uniprot_1"].isin(alzheimer_ids)
    & protein_info["uniprot_2"].isin(alzheimer_ids)
]

subgraph = nx.from_pandas_edgelist(
    alzheimer_edges,
    source="uniprot_1",
    target="uniprot_2",
    edge_attr="combined_score",
)

subgraph.add_nodes_from(alzheimer_ids)

# Gráfico estático de los 20 genes Alzheimer con mayor degree centrality
top_20_genes = gene_summary.sort_values(
    "degree_centrality", ascending=False
).head(20)
top_20_genes = top_20_genes.sort_values("degree_centrality")

plt.figure(figsize=(10, 8))
plt.barh(
    top_20_genes["symbol"],
    top_20_genes["degree_centrality"],
    color="#3182bd",
)
plt.xlabel("Degree centrality")
plt.ylabel("Gene")
plt.title("Top 20 genes de Alzheimer por degree centrality")
plt.tight_layout()
plt.savefig(FIGURE_PATH, dpi=300)
plt.close()

# Visualización interactiva con PyVis
net = Network(
    height="750px",
    width="100%",
    directed=False,
    notebook=False,
    bgcolor="#ffffff",
    font_color="#222222",
)

category_labels = {
    "hub": "Hub top 50 de G_giant",
    "chembl_target": "Alzheimer con diana ChEMBL",
    "no_chembl_target": "Alzheimer sin diana ChEMBL",
}

for index, row in gene_summary.iterrows():
    uniprot_id = row["uniprot_id"]
    category = row["category"]
    degree = round(row["degree_centrality"] * (len(giant_nodes) - 1))
    if row["chembl_target"]:
        target_html = escape(row["chembl_target"]).replace("; ", "<br>")
    else:
        target_html = "Ninguna"
    tooltip = (
        f"Gene symbol: {escape(row['symbol'])}"
        f"UniProt: {escape(uniprot_id)}"
        f"Global score: {row['globalScore']:.4f}"
        f"Degree en G_giant: {degree}"
        f"Registro ChEMBL de diana: "
        f"{'Sí' if row['chembl_target'] else 'No'}"
        f"Diana(s) ChEMBL: {target_html}<br>"
        f"Tipo: {category_labels[category]}"
    )
    net.add_node(
        uniprot_id,
        label=row["symbol"],
        title=tooltip,
        color=COLORS[category],
        value=max(degree, 1),
    )

# Añadir las interacciones al grafo de PyVis
for source, target, edge_data in subgraph.edges(data=True):
    net.add_edge(
        source,
        target,
        title=f"Combined score: {edge_data['combined_score']}",
    )

# Guardado del grafo interactivo
os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
net.write_html(OUTPUT_PATH, open_browser=False, notebook=False)

hub_count = gene_summary[gene_summary["category"] == "hub"].shape[0]

chembl_count = gene_summary[
    gene_summary["category"] == "chembl_target"
].shape[0]

no_chembl_count = gene_summary[
    gene_summary["category"] == "no_chembl_target"
].shape[0]

legend_items = (
    '<div style="margin:6px 0">'
    f'<span style="display:inline-block;width:12px;height:12px;'
    f'border-radius:50%;background:{COLORS["hub"]};'
    f'margin-right:7px"></span>{category_labels["hub"]} ({hub_count})</div>'
    '<div style="margin:6px 0">'
    f'<span style="display:inline-block;width:12px;height:12px;'
    f'border-radius:50%;background:{COLORS["chembl_target"]};'
    f'margin-right:7px"></span>{category_labels["chembl_target"]} '
    f"({chembl_count})</div>"
    '<div style="margin:6px 0">'
    f'<span style="display:inline-block;width:12px;height:12px;'
    f'border-radius:50%;background:{COLORS["no_chembl_target"]};'
    f'margin-right:7px"></span>{category_labels["no_chembl_target"]} '
    f"({no_chembl_count})</div>"
)

with open(OUTPUT_PATH, encoding="utf-8") as html_file:
    html_content = html_file.read()
html_content = html_content.replace(
    "</body>",
    (
        '<div style="position:fixed;right:16px;top:16px;z-index:1000;'
        'background:#fff;border:1px solid #bbb;border-radius:6px;'
        'padding:10px 14px;font:14px Arial,sans-serif">'
        "<b>Tipos de nodo</b>"
        f"{legend_items}</div></body>"
    ),
)

with open(OUTPUT_PATH, "w", encoding="utf-8") as html_file:
    html_file.write(html_content)

print(f"Genes Alzheimer incluidos: {len(subgraph.nodes)}")
print(f"Interacciones en el subgrafo: {len(subgraph.edges)}")
print(f"Grafo interactivo guardado en: {OUTPUT_PATH}")
print(f"Gráfico estático guardado en: {FIGURE_PATH}")

