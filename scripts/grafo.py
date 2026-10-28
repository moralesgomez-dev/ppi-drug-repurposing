import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os


PROJECT_ROOT = os.getcwd()
DATA_PATH = os.path.join(PROJECT_ROOT, "data")  
protein_info = pd.read_csv(os.path.join(DATA_PATH, "processed\\protein_info_csv"))
                            
protein_graph = nx.from_pandas_edgelist(protein_info,
                                        source="uniprot_1",
                                        target="uniprot_2",
                                        edge_attr="combined_score")
print("Nodos del grafo:", len(protein_graph.nodes))
print("Aristas del grafo:", len(protein_graph.edges))
print(nx.is_connected(protein_graph)) 
print(nx.number_connected_components(protein_graph))
# Tenemos 71 componentes diferentes en el grafo. Es decir, 71 bloques de componentes o proteinas que no se conectan entre si
# Tenemos que trabajar con el componente mas grande y sobre el calculas las diferentes metricas de centralidad

giant = max(nx.connected_components(protein_graph), key=len)
G_giant = protein_graph.subgraph(giant)
print(nx.is_connected(G_giant))
print(nx.number_connected_components(G_giant))
# Como vemos ya nos hemos quedado solo con un componente
print("Nodos del grafo:", len(G_giant.nodes))
print("Aristas del grafo:", len(G_giant.edges))
# Este componente principal resulta ser casi la red entera 

FIGURES_PATH = os.path.join(PROJECT_ROOT, "results", "figures")
os.makedirs(FIGURES_PATH, exist_ok=True)

# Metricas de centralidad:
degree_cent = nx.degree_centrality(G_giant)
# betweenness_cent = nx.betweenness_centrality(G_giant, k=500)
# closeness_cent = nx.closeness_centrality(G_giant, k=500)

# Nos quedamos solo con los top 50 hubs por degree
top_50_hubs = sorted(degree_cent, key=degree_cent.get, reverse=True)[:50]
print("Top 50 hubs por degree:", top_50_hubs)

# Distribucion de grados
G_giant_degrees = [d for n, d in G_giant.degree()]
plt.figure(figsize=(10, 6))
plt.hist(G_giant_degrees, bins=50)
plt.xlabel("Grado")
plt.ylabel("Frecuencia")
plt.title("Distribución de Grados")
plt.tight_layout()
plt.savefig(
    os.path.join(FIGURES_PATH, "degree_distribution.png"),
    dpi=300,
)
plt.close()

# Distribucion de centralidad de grado
degree_centrality_values = list(degree_cent.values())
plt.figure(figsize=(10, 6))
plt.hist(degree_centrality_values, bins=50)
plt.xlabel("Centralidad de Grado")
plt.ylabel("Frecuencia")
plt.title("Distribución de Centralidad de Grado")
plt.tight_layout()
plt.savefig(
    os.path.join(FIGURES_PATH, "degree_centrality_distribution.png"),
    dpi=300,
)
plt.close()

# Subgrafo de los top 50 hubs con Pyvis - HTML interactivo
from pyvis.network import Network

top_50_subgraph = G_giant.subgraph(top_50_hubs).copy()
network = Network(
    height="750px",
    width="100%",
    directed=False,
    notebook=False,
    bgcolor="#ffffff",
    font_color="#222222",
)
network.from_nx(top_50_subgraph)

for node in network.nodes:
    protein = node["id"]
    node["label"] = str(protein)
    node["title"] = f"Proteína: {protein} - Grado: {top_50_subgraph.degree[protein]}"

for edge in network.edges:
    if "combined_score" in edge:
        edge["title"] = f"Puntuación combinada: {edge['combined_score']}"

network.show_buttons(filter_=["physics"])
output_path = os.path.join(FIGURES_PATH, "top_50_hubs.html")
network.write_html(output_path, open_browser=False, notebook=False)
print(f"Grafo interactivo guardado en: {output_path}")

degree_df = pd.DataFrame(
    degree_cent.items(),
    columns=["uniprot_id", "degree_centrality"]
)

degree_df.head()
degree_df.sample(5)

degree_df.to_csv(os.path.join(DATA_PATH, "processed\\degree_centrality.csv"), index=False)

