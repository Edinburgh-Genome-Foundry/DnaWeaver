from bokeh.io import output_notebook, output_file
from bokeh.plotting import figure, show, ColumnDataSource, hplot
from bokeh.models import FixedTicker, Range1d, TapTool, OpenURL, Rect, CustomJS, HoverTool
import pandas as pd
import matplotlib.cm as cm
import numpy as np

from bokeh.io import output_notebook, output_file
from bokeh.plotting import figure, show, ColumnDataSource, hplot
from bokeh.models import FixedTicker, Range1d, TapTool, OpenURL, Rect, CustomJS, HoverTool
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches


def plot_assembly_tree(tree, color_palette=None, parts_offset=0,
                       backend="bokeh", plot_top_assembly=True,
                       ax=None):
    rectangles = []

    def rec(node, depth=0, xstart=0, xend=None):
        quote, children = node
        if plot_top_assembly or (depth > 0):
            rectangles.append({
                "top": -depth + .5 + (not plot_top_assembly),
                "bottom": -depth + (not plot_top_assembly),
                "left": xstart,
                "right": xend,
                "source": str(quote.source),
                "lead_time": str(quote.lead_time),
                "price": str(quote.price),
                "has_children": children != {}
            })
        for (start, end), child in children.items():
            rec(child, depth=depth + 1, xstart=xstart + start,
                xend=xstart + end)
    rec(tree, xend=len(tree[0].sequence))
    tops = list(set([r["top"] for r in rectangles]))
    offsets = {top: 0 for top in tops}

    if color_palette is None:
        sources = list(set([r["source"] for r in rectangles]))
        L = len(sources)
        colors_cycle = [
            '#%02x%02x%02x' %
            tuple(
                (255 * np.array(cm.Paired(0.21 * i % 1.0)[:3])).astype(int))
            for i in range(L)]
        color_palette = {source: color
                         for source, color in zip(sources, colors_cycle)}

    guides = []
    for top in sorted(tops):
        rects = sorted([r for r in rectangles if r["top"] == top],
                       key=lambda r: r["left"])
        for i, rect in enumerate(rects):
            color = ["#e5ecff", "#ffffe5"][i % 2]
            if rect["has_children"]:
                guides.append({
                    "top": rect["top"],
                    "bottom": rect["top"] - 1,
                    "left": rect["left"],
                    "right": rect["right"],
                    "color": color
                })

    for rectangle in sorted(rectangles, key=lambda r: r["left"]):
        rectangle["color"] = color_palette[rectangle["source"]]
        if parts_offset:
            offset = offsets[rectangle["top"]]
            offsets[rectangle["top"]] = (offsets[rectangle["top"]] + 1) % 2
            rectangle["top"] += parts_offset * offset
            rectangle["bottom"] += parts_offset * offset

    if backend == "bokeh":
        return _bokeh_plot(rectangles, guides, tops)
    elif backend == "matplotlib":
        return _matplotlib_plot(rectangles, guides, tops, ax)
    else:
        raise ValueError("Backend should be bokeh or matplotlib")


def _bokeh_plot(rectangles, guides, tops):
    hover = HoverTool(
        always_active=True,
        tooltips="""
        <div> <b> @source </b><br/>
        @price $ @lead_time days
        </div>
    """,
        names=["parts"]
    )

    L = max([r["right"] for r in rectangles])
    p = figure(tools=[hover, "xwheel_zoom,xpan,reset"],
               responsive=True,
               plot_height=300,
               x_range=Range1d(int(-0.05 * L), int(1.05 * L)),
               y_range=Range1d(min(tops) - 1, -min(tops)),
               logo=None)
    source_guide = ColumnDataSource(pd.DataFrame.from_records(guides))
    p.quad(name="guides",
           top="top", bottom="bottom", left="left",
           right="right", color="color", source=source_guide, alpha=0.5)
    dataframe = pd.DataFrame.from_records(rectangles)
    for sourcename, df in dataframe.groupby("source"):
        source = ColumnDataSource(df)
        print sourcename, len(df)
        p.quad(name="parts",
               top="top", bottom="bottom", left="left",
                   right="right", color="color", source=source, line_width=1,
                   line_color="black", alpha=1, legend=sourcename)

    p.yaxis.visible = None
    p.outline_line_color = None
    p.grid.grid_line_color = None
    return p


def _matplotlib_plot(rectangles, guides, tops, ax):
    L = max([r["right"] for r in rectangles])
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(8, 8))
    ax.set_xlim((-0.05 * L, 1.05 * L))
    ax.set_ylim((min(tops) - 1, 1.5))

    for g in guides:
        patch = patches.Rectangle((g["left"], g["bottom"]),
                                  g["right"] - g["left"],
                                  g["top"] - g["bottom"],
                                  color=g["color"],
                                  ec="none")
        ax.add_patch(patch)

    seen_sources = set([])
    legend_handles = []
    for g in rectangles:
        if g["source"] not in seen_sources:
            legend_patch = patches.Patch(facecolor=g["color"],
                                         label=g["source"],
                                         edgecolor="k",
                                         linewidth=1)
            legend_handles.append(legend_patch)
            seen_sources.add(g["source"])
            plt.legend(handles=[legend_patch])
        patch = patches.Rectangle((g["left"], g["bottom"]),
                                  g["right"] - g["left"],
                                  g["top"] - g["bottom"],
                                  facecolor=g["color"],
                                  edgecolor="k", linewidth=1)
        ax.add_patch(patch)
    ax.axis("off")
    return ax, legend_handles
