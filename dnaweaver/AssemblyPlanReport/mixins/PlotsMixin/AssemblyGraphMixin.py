from collections import defaultdict
from copy import deepcopy
import matplotlib.font_manager as fm
import numpy as np
from ...config import SETTINGS
from .plot_tree_graph import plot_tree_graph


class AssemblyGraphMixin:
    def plot_assembly_graph(self, ax=None, margin=None, textprops=None, scale=1.0):
        """Plot the complete assembly graph.

        Returns
        -------
        elements_positions, ax
        Dictionary of element positions, matplotlib ax.
        """

        nodes_dict = {}
        levels = defaultdict(lambda *a: [])
        edges = []
        tree = deepcopy(self.plan)

        def rec(node, depth=0):
            if node.get("_visited", False):
                return
            nodes_dict[node.id] = node
            node["_visited"] = True
            assembly_plan = node.pop("assembly_plan")
            levels[depth].append(node.id)
            for other in assembly_plan:
                edges.append([other.id, node.id])
                rec(other, depth + 1)

        rec(tree)
        levels = [levels[i] for i in range(max(levels) + 1)][::-1]

        fontawesome = fm.FontProperties(
            fname=SETTINGS["fontawesome-ttf-path"],
            size=13 * scale,
            family="sans-serif",
        )
        if textprops is None:
            textprops = fm.FontProperties(
                fname=SETTINGS["OpenSans-ttf-path"],
                size=12 * scale,
                family="sans-serif",
            )

        def draw_node(x, y, node_id, ax):
            node = nodes_dict[node_id]
            icon = self.sources[node.source]._report_fa_symbol
            ax.text(
                x,
                y,
                node_id,
                horizontalalignment="left",
                verticalalignment="center",
                fontproperties=textprops,
            )
            ax.text(
                x - 0.01 * np.sqrt(scale),
                y,
                icon,
                horizontalalignment="right",
                verticalalignment="center",
                fontproperties=fontawesome,
            )

        all_elements = sorted(sum(levels, []))
        ypos = {
            el: 1.0 * (i + 1) / (len(all_elements) + 2)
            for i, el in enumerate(all_elements)
        }
        for el in all_elements:
            children = [e2 for (e2, e1) in edges if e1 == el]
            if children != []:
                ypos[el] = 1.0 * sum(ypos[e] for e in children) / len(children)

        xpos = {
            el: 1.0 * (1 + x) / (len(levels) + 1)
            for x, elements in enumerate(levels)
            for el in elements
        }

        elements_positions = {el: (xpos[el], ypos[el]) for el in all_elements}
        return plot_tree_graph(
            levels,
            edges,
            draw_node,
            elements_positions=elements_positions,
            ax=ax,
            edge_left_space=0.06,
            edge_right_space=0.03,
            margin=margin,
            height_factor=0.40,
            width_factor=5.5,
            scale=scale,
        )
