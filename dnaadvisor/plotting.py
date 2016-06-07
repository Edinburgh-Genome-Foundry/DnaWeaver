from bokeh.io import output_notebook, output_file
from bokeh.plotting import figure, show, ColumnDataSource
from bokeh.models import FixedTicker, Range1d, TapTool, OpenURL, Rect, CustomJS, HoverTool
import pandas as pd
import matplotlib.cm as cm
import numpy as np


def plot_ordering_tree(tree, color_palette=None):
    rectangles = []
        def rec(node, depth=0, xstart=0, xend=None):
            quote, children = node
            rectangles.append({
                "top": -depth + .5,
                "bottom": -depth,
                "left": xstart,
                "right": xend,
                "source": str(quote.source),
                "lead time": str(quote.lead_time),
                "price": str(quote.price),
                "has_children": children != {}
            })
            for (start, end), child in children.items():
                rec(child, depth=depth + 1, xstart=xstart +
                    start, xend=xstart + end)
        rec(tree)
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
        source_guide = ColumnDataSource(pd.DataFrame.from_records(guides))

        for rectangle in sorted(rectangles, key=lambda r: r["left"]):
            rectangle["color"] = color_palette[rectangle["source"]]
            if parts_offset:
                offset = offsets[rectangle["top"]]
                offsets[rectangle["top"]] = (offsets[rectangle["top"]] + 1) % 2
                rectangle["top"] += parts_offset * offset
                rectangle["bottom"] += parts_offset * offset

        source = ColumnDataSource(pd.DataFrame.from_records(rectangles))

        hover = HoverTool(
            always_active=True,
            tooltips="""
            <div> <b> @source </b><br/>
            @price $ @lead_time days
            </div>
        """,
            names=["parts"]
        )
        responsive = True
        L = len(sequence)
        p = figure(tools=[hover, "xwheel_zoom,xpan,reset"],
                   responsive=True,
                   plot_height=100,
                   x_range=Range1d(int(-0.05 * L), int(1.05 * L)),
                   y_range=Range1d(min(tops) - 1, 0),
                   logo=None)

        p.quad(name="guides",
               top="top", bottom="bottom", left="left",
               right="right", color="color", source=source_guide, alpha=0.5)

        p.quad(name="parts",
               top="top", bottom="bottom", left="left",
                   right="right", color="color", source=source, line_width=1,
                   line_color="black", alpha=1)

        p.yaxis.visible = None
        p.outline_line_color = None
        p.grid.grid_line_color = None
        return p
