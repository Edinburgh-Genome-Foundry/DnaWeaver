from dna_features_viewer import GraphicRecord
import matplotlib.pyplot as plt
from copy import deepcopy
import matplotlib.patches as mpatches
import matplotlib.font_manager as fm
from Bio import SeqIO
from ...config import SETTINGS


class AssemblyBlocksMixin:
    def plot_assembly_blocks(
        self,
        parts_offset=0,
        plot_top_assembly=True,
        ax=None,
        edge_widths=None,
        legend=False,
        legend_offset=-0.05,
    ):
        """Return a Matplotlib or Bokeh plot of the assembly tree of blocks.

        Parameters
        ----------

        parts_offset
          Offset applied so that consecutive blocks are not exactly on the same
          level. Can go from 0 (flat line of blocks) to e.g. 1.

        plot_top_assembly
          Whether the top assembly (which is just one big block) should be
          plotted or not.

        ax
          A Matplotlib Axes object. If no ax is provided, a new figure and ax
          are generated.

        edge_widths
          A dict {sourcename : width} indicating the widths of the rectangles
          based on the source name. Rectangles with a very small width will be
          edgeless.

        legend
          Whether the legend is included in the ax.
        """
        rectangles = []
        if edge_widths is None:
            edge_widths = {}
        if not hasattr(self, "plan"):
            self.compute_full_assembly_plan()
        tree = deepcopy(self.plan)

        def rec(_quote, depth=0, xstart=0, xend=None):
            children = _quote.pop("assembly_plan")
            # print _quote.final_location
            if _quote.final_location is not None:
                left, right = _quote.final_location
            else:
                left, right = xstart, xend

            if children is None:
                children = []
            if plot_top_assembly or (depth > 0):
                source = self.sources[_quote.source]
                rectangles.append(
                    {
                        "top": -depth + 0.5 + (not plot_top_assembly),
                        "bottom": -depth + (not plot_top_assembly),
                        "left": left,
                        "right": right,
                        "source": str(_quote.source),
                        "lead_time": str(_quote.lead_time),
                        "price": str(_quote.price),
                        "length": len(_quote.sequence),
                        "has_children": children != [],
                        "line_width": 1,
                        "operation_type": source.operation_type,
                        "fa_symbol": source._report_fa_symbol_plain,
                        "color": source._report_color,
                    }
                )
            for child in children:
                rec(
                    child,
                    depth=depth + 1,
                    xstart=xstart + child.segment_start,
                    xend=xstart + child.segment_end,
                )

        rec(tree, xend=len(tree.sequence))
        tops = list(set([r["top"] for r in rectangles]))
        offsets = {top: 0 for top in tops}

        guides = []
        for top in sorted(tops):
            rects = sorted(
                [r for r in rectangles if r["top"] == top], key=lambda r: r["left"]
            )
            for i, rect in enumerate(rects):
                color = ["#e5ecff", "#ffffe5"][i % 2]
                if rect["has_children"]:
                    guides.append(
                        {
                            "top": rect["top"],
                            "bottom": rect["top"] - 1,
                            "left": rect["left"],
                            "right": rect["right"],
                            "color": color,
                        }
                    )
        for rectangle in sorted(rectangles, key=lambda r: r["left"]):
            if parts_offset:
                offset = offsets[rectangle["top"]]
                offsets[rectangle["top"]] = (offsets[rectangle["top"]] + 1) % 2
                rectangle["top"] += parts_offset * offset
                rectangle["bottom"] += parts_offset * offset

        if hasattr(self, "genbank") and (self.genbank is not None):
            record = SeqIO.read(self.genbank, "genbank")
        else:
            record = None
        return _matplotlib_plot_assembly_blocks(
            rectangles,
            guides,
            tops,
            ax,
            record=record,
            legend=legend,
            legend_offset=legend_offset,
        )


def _matplotlib_plot_assembly_blocks(
    rectangles,
    guides,
    tops,
    ax=None,
    fig_height="auto",
    fig_width=8,
    legend=False,
    textprops=None,
    record=None,
    legend_offset=-0.05,
):
    """Plot the assembly block rectangles using matplotlib."""

    if record is not None:
        fig, record_ax, blocks_ax = plt.subplots(
            2, 1, figsize=(7, 7), gridspec_kw={"height_ratios": [4, 1]}
        )
        grecord = GraphicRecord.from_biopython_record(
            record, fun_color=lambda a: "#aabbff"
        )
        grecord.plot(record_ax)
        res = _matplotlib_plot_assembly_blocks(
            rectangles,
            guides,
            tops,
            ax=blocks_ax,
            fig_height="auto",
            fig_width=8,
            legend=False,
            textprops=None,
            record=None,
        )
        record_ax.set_xlim(blocks_ax.get_xlim())
        return res

    L = max([r["right"] for r in rectangles])

    if textprops is None:
        textprops = fm.FontProperties(
            fname=SETTINGS["OpenSans-ttf-path"], size=12, family="sans-serif"
        )

    if ax is None:
        if fig_height == "auto":
            fig_height = 8
        fig, ax = plt.subplots(1, figsize=(8, 8))
    ax.set_xlim((-0.05 * L, 1.05 * L))
    ax.set_ylim((min(tops) - 1, 1.5))

    for g in guides:
        patch = mpatches.Rectangle(
            (g["left"], g["bottom"]),
            g["right"] - g["left"],
            g["top"] - g["bottom"],
            color=g["color"],
            ec="none",
        )
        ax.add_patch(patch)

    seen_sources = set([])
    legend_handles = []
    sorted_rectangles = sorted(rectangles, key=lambda r: (-r["bottom"], r["left"]))
    for g in sorted_rectangles:
        if g["source"] not in seen_sources:
            legend_patch = mpatches.Patch(
                facecolor=g["color"], label=g["source"], edgecolor="k", linewidth=1.0,
            )
            legend_handles.append(legend_patch)
            seen_sources.add(g["source"])
            if legend:
                plt.legend(handles=[legend_patch])
        # linewidth = g["line_width"]
        width = g["right"] - g["left"]
        line_width = 1.0 if (1.0 * width / L) > 0.002 else 0
        patch = mpatches.Rectangle(
            (g["left"], g["bottom"]),
            width,
            g["top"] - g["bottom"],
            facecolor=g["color"],
            edgecolor="k",
            linewidth=line_width,
        )
        ax.add_patch(patch)
    ax.axis("off")

    if legend:
        if ax.legend_ is not None:
            ax.legend_.remove()
        legend = ax.legend(
            handles=legend_handles,
            frameon=False,
            ncol=2,
            loc=2,
            bbox_to_anchor=(0.0, legend_offset),
        )
        ltext = legend.get_texts()
        plt.setp(ltext, fontproperties=textprops)
        ax.figure.subplots_adjust(hspace=0.0)

    return ax, legend_handles
