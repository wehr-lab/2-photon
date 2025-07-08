import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Functions to compute metrics
def compute_rich_club(G, k=5):
    rich_nodes = [n for n in G.nodes if G.degree[n] > k]
    subgraph = G.subgraph(rich_nodes)
    possible_edges = len(rich_nodes) * (len(rich_nodes) - 1) / 2
    return subgraph.number_of_edges() / possible_edges if possible_edges else 0

def compute_small_worldness(G):
    if not nx.is_connected(G):
        largest_cc = max(nx.connected_components(G), key=len)
        G = G.subgraph(largest_cc).copy()
    C = nx.average_clustering(G)
    try:
        L = nx.average_shortest_path_length(G)
    except nx.NetworkXError:
        return 0
    deg_seq = [d for _, d in G.degree()]
    G_rand = nx.configuration_model(deg_seq)
    G_rand = nx.Graph(G_rand)
    G_rand.remove_edges_from(nx.selfloop_edges(G_rand))
    if not nx.is_connected(G_rand):
        largest_cc = max(nx.connected_components(G_rand), key=len)
        G_rand = G_rand.subgraph(largest_cc).copy()
    C_rand = nx.average_clustering(G_rand)
    try:
        L_rand = nx.average_shortest_path_length(G_rand)
    except nx.NetworkXError:
        return 0
    return (C / C_rand) / (L / L_rand) if C_rand and L_rand else 0

# Simulation settings
n_nodes = 500
n_subjects = 5
attack_steps = 15

# Store metric trajectories
rc_rich, rc_rand = [], []
sw_rich, sw_rand = [], []

bin_mat = np.load("/Users/praveslamichhane/Desktop/network_comp/data/bin_mat.npy", allow_pickle=True)
# Simulation loop
for _ in range(n_subjects):
    
    G_base = nx.from_numpy_array(bin_mat)

    G1 = G_base.copy()
    G2 = G_base.copy()
    r1, r2 = [], []
    s1, s2 = [], []

    for _ in range(attack_steps):
        r1.append(compute_rich_club(G1))
        s1.append(compute_small_worldness(G1))
        r2.append(compute_rich_club(G2))
        s2.append(compute_small_worldness(G2))
        G1.remove_node(max(G1.degree, key=lambda x: x[1])[0])
        G2.remove_node(np.random.choice(list(G2.nodes)))

    rc_rich.append(r1)
    rc_rand.append(r2)
    sw_rich.append(s1)
    sw_rand.append(s2)

# Prepare DataFrames
def aggregate_results(trajectories, label):
    df = pd.DataFrame(trajectories).T
    df['Step'] = df.index
    df = df.melt(id_vars='Step', var_name='Subject', value_name='Value')
    df['AttackType'] = label
    return df

df_rc = pd.concat([
    aggregate_results(rc_rich, 'Rich-Club'),
    aggregate_results(rc_rand, 'Random')
])
df_rc['Metric'] = 'Rich-Club Coefficient'

df_sw = pd.concat([
    aggregate_results(sw_rich, 'Rich-Club'),
    aggregate_results(sw_rand, 'Random')
])
df_sw['Metric'] = 'Small-Worldness'

df_all = pd.concat([df_rc, df_sw], ignore_index=True)

# Plotting
sns.set(style="whitegrid", font_scale=1.2)
metrics = df_all['Metric'].unique()
fig, axes = plt.subplots(1, len(metrics), figsize=(14, 6), sharex=True)

for ax, metric in zip(axes, metrics):
    sns.lineplot(
        data=df_all[df_all['Metric'] == metric],
        x='Step', y='Value',
        hue='AttackType', ci=95,
        palette='Set1', linewidth=2.5, ax=ax
    )
    ax.set_title(metric)
    ax.set_xlabel("Nodes Removed")
    ax.set_ylabel("Value")
    ax.legend(title='Attack Type')

plt.suptitle("Attack Impact on Network Metrics (95% Confidence Interval)", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()
