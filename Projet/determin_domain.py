import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.widgets import Button

# Chargement des données
file_path = "data/mesh.txt"  # Remplace par le chemin de ton fichier
with open(file_path, "r") as f:
    lines = f.readlines()

# Extraction des nœuds
nodes = {}
n_edges = []
triangles = []
domains = {}
current_section = None

domain_ids = []
current_domain_index = 0
current_domain = None
parsing_elements = False

for line in lines:
    if "Number of nodes" in line:
        current_section = "nodes"
    elif "Number of edges" in line:
        current_section = "edges"
    elif "Number of tiangles" in line:
        current_section = "triangles"
    elif "Domain" in line:
        current_section = "domains"
        parsing_elements = False
        domain_id = int(line.split(":")[1].strip())
        domains[domain_id] = []
        domain_ids.append(domain_id)
        current_domain = domain_id
    elif "Number of elements" in line:
        parsing_elements = True
    elif parsing_elements and line.strip():
        elements = list(map(int, line.split()))
        domains[current_domain].extend(elements)
    elif current_section == "nodes" and ":" in line:
        parts = line.split(":")
        idx = int(parts[0].strip())
        x, y = map(float, parts[1].split())
        nodes[idx] = (x, y)
    elif current_section == "edges" and ":" in line:
        parts = line.split(":")
        edge_data = parts[1].split()
        if len(edge_data) == 2:
            n1, n2 = map(int, edge_data)
            n_edges.append((n1, n2))
    elif current_section == "triangles" and ":" in line:
        parts = line.split(":")
        n1, n2, n3  = map(int, parts[1].split())
        triangles.append((n1, n2, n3))

# Création du graphe
print(f"Nombre de nœuds: {len(nodes)}")
print(f"Nombre d'arêtes: {len(n_edges)}")
print(f"Nombre de triangles: {len(triangles)}")
print(f"Nombre de domaines: {len(domains)}")
# Création du graphe
G = nx.Graph()
G.add_nodes_from(nodes.keys())
G.add_edges_from(n_edges)

# Configuration de la figure
fig, ax = plt.subplots(figsize=(10, 10))
pos = {k: (v[0], -v[1]) for k, v in nodes.items()}  # Inverser l'axe Y pour correspondre au repère usuel

def update_plot():
    ax.clear()
    
    current_domain = domain_ids[current_domain_index]
    print(current_domain)
    domain_nodes = domains[current_domain]
    domain_edges = [(n1, n2) for n1, n2 in n_edges if n1 in domain_nodes and n2 in domain_nodes]

    # Vérification de la validité des nœuds dans le domaine
    invalid_nodes = [n for n in domain_nodes if n not in nodes or n == -1]
    if invalid_nodes:
        print(f"Avertissement : Nœuds invalides dans le domaine {current_domain}: {invalid_nodes}")

    # Filtrer les nœuds valides
    valid_nodes = [n for n in domain_nodes if n in nodes and n != -1]
    valid_pos = {n: pos[n] for n in valid_nodes}
    
    # Tracer toutes les arêtes en gris
    nx.draw(G, pos, node_size=1, edge_color='lightgray', alpha=0.5, ax=ax, with_labels=False)
    
    # Tracer les arêtes du domaine en rouge
    H = nx.Graph()
    H.add_edges_from(domain_edges)
    nx.draw(H, pos, node_size=5, edge_color='red', alpha=0.9, ax=ax, with_labels=False)
    
    # Mettre en évidence les nœuds du domaine
    x = [valid_pos[e][0] for e in valid_nodes]
    y = [valid_pos[e][1] for e in valid_nodes]
    ax.scatter(x, y, color='red', s=30, label=f'Domaine {current_domain}', alpha=0.9)
    
    ax.set_title(f"Visualisation du Domaine {current_domain}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.legend()
    plt.draw()



def next_domain(event):
    global current_domain_index
    current_domain_index = (current_domain_index + 1) % len(domain_ids)
    update_plot()

def prev_domain(event):
    global current_domain_index
    current_domain_index = (current_domain_index - 1) % len(domain_ids)
    update_plot()

# Boutons pour changer de domaine
axprev = plt.axes([0.7, 0.01, 0.1, 0.05])
axnext = plt.axes([0.81, 0.01, 0.1, 0.05])
btn_prev = Button(axprev, "<-")
btn_next = Button(axnext, "->")
btn_prev.on_clicked(prev_domain)
btn_next.on_clicked(next_domain)

update_plot()
plt.show()