"""Nice function to make figures for slides (not used in the framework)."""

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def plot_decomposition_graph(
    graph,
    nodes_color="#6886b7",
    weight="weight",
    colormap="jet",
    ax=None,
    edge_width=1,
    edges_step=1,
    edge_alpha=0.2,
    figsize=(8, 8),
):
    xx = np.array(sorted(graph.nodes))
    L = xx.max()
    edges = [
        (start, end, data[weight] / (end - start))
        for start, end, data in graph.edges(data=True)
    ]
    edges = edges[::edges_step]
    edges = sorted(edges, key=lambda e: e[2])
    max_segment_length = max([end - start for (start, end, _) in edges])
    weights = np.array([w for (_, _, w) in edges])
    normalized_weights = (255 * weights / weights.max()).astype("uint8")
    colormap = cm.__dict__[colormap]
    colors = colormap(normalized_weights, alpha=edge_alpha)
    if ax is None:
        fig, ax = plt.subplots(1, figsize=figsize)
    ax.axis("off")
    #     print (list(colors))
    for (start, end, w), color in zip(edges, colors):
        xc = 0.5 * (start + end)
        half = 0.5 * abs(end - start)
        ax.add_patch(
            mpatches.Arc(
                (xc, 0),
                2 * half,
                2 * half,
                theta1=0,
                theta2=180,
                facecolor="none",
                ls="-",
                edgecolor=color,
                linewidth=1,
            )
        )
    ax.plot(xx, [0 for x in xx], marker="o", c=nodes_color)
    ax.set_aspect("equal")
    ax.set_xlim(-1, L + 1)
    ax.set_ylim(-1, max_segment_length / 2 + 1)
