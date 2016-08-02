#!/usr/bin/env python
# -*- coding: utf-8 -*-

import textwrap

from bokeh.io import output_notebook, output_file
from bokeh.plotting import figure, show, ColumnDataSource, hplot
from bokeh.models import (FixedTicker, Range1d, TapTool, OpenURL, Rect,
                          CustomJS, HoverTool)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.font_manager as fm
from matplotlib.path import Path

import config


def plot_assembly_tree(tree, color_palette=None, parts_offset=0,
                       backend="bokeh", plot_top_assembly=True,
                       ax=None, edge_widths=None, legend=False):
    rectangles = []
    if edge_widths is None:
        edge_widths = {}

    def rec(node, depth=0, xstart=0, xend=None):
        quote, children = node.quote, node.segments
        if plot_top_assembly or (depth > 0):
            rectangles.append({
                "top": -depth + .5 + (not plot_top_assembly),
                "bottom": -depth + (not plot_top_assembly),
                "left": xstart,
                "right": xend,
                "source": str(quote.source),
                "lead_time": str(quote.lead_time),
                "price": str(quote.price),
                "has_children": children != {},
                "line_width": edge_widths.get(str(quote.source), 1)
            })
        for (start, end), child in children.items():
            rec(child, depth=depth + 1, xstart=xstart + start,
                xend=xstart + end)
    rec(tree, xend=len(tree.quote.sequence))
    tops = list(set([r["top"] for r in rectangles]))
    offsets = {top: 0 for top in tops}

    if color_palette is None:
        sources = list(set([r["source"] for r in rectangles]))
        L = len(sources)
        colors_cycle = [
            '#%02x%02x%02x' %
            tuple((255*(0.4 +
                      0.6*np.array(cm.Paired(0.21 * i % 1.0)[:3]))).astype(int))
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
        return _bokeh_plot_assembly_tree(rectangles, guides, tops)
    elif backend == "matplotlib":
        return _matplotlib_plot_assembly_tree(rectangles, guides, tops, ax,
                                              legend=legend)
    else:
        raise ValueError("Backend should be bokeh or matplotlib")


def _bokeh_plot_assembly_tree(rectangles, guides, tops):
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
        p.quad(name="parts",
               top="top", bottom="bottom", left="left",
                   right="right", color="color", source=source, line_width=1,
                   line_color="black", alpha=1, legend=sourcename)

    p.yaxis.visible = None
    p.outline_line_color = None
    p.grid.grid_line_color = None
    return p


def _matplotlib_plot_assembly_tree(rectangles, guides, tops, ax, legend=False,
                                   textprops=None):
    L = max([r["right"] for r in rectangles])

    if textprops is None:
        textprops = fm.FontProperties(
            fname=config.SETTINGS["OpenSans-ttf-path"],
            size=12, family="sans-serif"
        )

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
                                         linewidth=1.0)
            legend_handles.append(legend_patch)
            seen_sources.add(g["source"])
            plt.legend(handles=[legend_patch])
        # linewidth = g["line_width"]
        width = g["right"] - g["left"]
        line_width = 1.0 if (1.0*width/L)>.01 else 0
        patch = patches.Rectangle((g["left"], g["bottom"]),
                                  width,
                                  g["top"] - g["bottom"],
                                  facecolor=g["color"],
                                  edgecolor="k",
                                  linewidth=line_width)
        ax.add_patch(patch)
    ax.axis("off")

    if legend:
        ax.legend_.remove()
        legend = ax.legend(handles=legend_handles,
                           frameon=False, ncol=2, loc=2,
                           bbox_to_anchor=(0.0, -.05))
        ltext = legend.get_texts()
        plt.setp(ltext, fontproperties=textprops)
        ax.figure.subplots_adjust(hspace=.0)

    return ax, legend_handles


def plot_supply_graph(top_source, ax=None, icons="default", textprops=None,
                      scale=1.0):

    edges, levels = top_source.compute_supply_graph()

    if icons == "default":
        icons = {
            "PcrOutStation": u"",
            "ExternalDnaOffer": u"",
            "DnaAssemblyStation": u"",
            "DnaSourcesComparator": u"",
            "PartsLibrary": u"",
        }

    fontawesome = fm.FontProperties(
        fname=config.SETTINGS["fontawesome-ttf-path"],
        size=18*scale, family="sans-serif"
    )

    if textprops is None:
        textprops = fm.FontProperties(
            fname=config.SETTINGS["OpenSans-ttf-path"],
            size=12*scale, family="sans-serif"
        )

    source_positions = {}

    if ax is None:
        width = 2.5 * len(levels)*scale
        height = 2 * max([len(lvl) for lvl in levels])*scale
        fig, ax = plt.subplots(1, figsize=(width, height))

    for lvl, sources in enumerate(levels):
        yy = np.linspace(0, 1, len(sources) + 2)[1:-1]
        x = 1.0 * (1 + lvl) / (len(levels) + 1)
        for y, source in zip(yy, sources):
            source_positions[source] = (x, y)
            icon = icons[source.__class__.__name__]
            ax.text(x, y, icon, horizontalalignment="center",
                    verticalalignment="center", fontproperties=fontawesome)
            if hasattr(source, "name"):
                text = "\n".join(textwrap.wrap(source.name, 13))

                ax.text(x, y + 0.035*scale**.3, text,
                        horizontalalignment="center",
                        verticalalignment="bottom", fontproperties=textprops)

    for source1, source2 in edges:
        x1, y1 = source_positions[source1]
        x2, y2 = source_positions[source2]
        x1 += 0.015*scale
        x2 += -0.015*scale
        ax.add_patch(
            patches.PathPatch(
                Path([(x1, y1), (0.9 * x2 + 0.1 * x1, y1),
                      (0.1 * x2 + 0.9 * x1, y2), (x2, y2)],
                     [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                facecolor='none', lw=1*scale
            )
        )

    ax.axis("off")
    return ax


def _matplotlib_plot_assembly_timeline(assembly_tree, deadline=None, ax=None,
                                       rectangle_color="#bbbbff", scale=1.0):
    if deadline is None:
        deadline = assembly_tree.quote.lead_time

    assembly_tree.propagate_deadline(deadline)
    assemblies_list = assembly_tree.return_tree_as_list()
    assemblies_list = sorted(assemblies_list, key=lambda a: a.id)

    textprops = fm.FontProperties(
        fname=config.SETTINGS["OpenSans-ttf-path"],
        size=12 * scale, family="sans-serif"
    )
    fontawesome = fm.FontProperties(
        fname=config.SETTINGS["fontawesome-ttf-path"],
        size=12 * scale, family="sans-serif"
    )
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(16 * scale,
                                           len(assemblies_list) * 0.3 * scale))
    positions = {}
    icons = {
        "PcrOutStation": u"",
        "ExternalDnaOffer": u"",
        "DnaAssemblyStation": u"",
        "DnaSourcesComparator": u"",
        "PartsLibrary": u"",
    }

    if textprops is None:
        textprops = fm.FontProperties(
            fname=config.SETTINGS["OpenSans-ttf-path"],
            size=12, family="sans-serif"
        )

    for i, assembly in enumerate(assemblies_list):
        source = assembly.quote.source
        color = "k"
        if source.__class__.__name__ == "ExternalDnaOffer":
            color = "#888888"
        x2 = assembly.deadline
        x1 = assembly.deadline - assembly.step_duration
        positions[assembly.id] = ([x1, i], [x2, i])
        patch = patches.Rectangle(
            (assembly.deadline - assembly.step_duration, i - .2),
            assembly.step_duration, .4, color=rectangle_color,
            ec="none"
        )

        ax.plot([0, assembly.deadline], [i, i], c="k", alpha=0.3, lw=.5*scale,
                zorder=-1000, ls="-")
        ax.add_patch(patch)

        ax.text(-1.5, i, assembly.id, horizontalalignment="right",
                verticalalignment="center", fontproperties=textprops,
                color=color)
        icon = icons[source.__class__.__name__]
        ax.text(-.75, i, icon, horizontalalignment="center",
                verticalalignment="center", fontproperties=fontawesome,
                color=color)

    for assembly in assemblies_list:
        xa, ya = positions[assembly.id][0]
        for segment, child in assembly.segments.items():
            xc, yc = positions[child.id][1]
            ax.plot([xc, xa], [yc, ya], lw=1 * scale, color=color)

    ax.set_ylim((-1, len(assemblies_list) + 1))
    ax.set_xlim(-2, deadline)
    ax.set_frame_on(False)
    ax.yaxis.set_visible(False)
    ticks = list(np.arange(0, deadline, deadline / 6).astype(int))

    ax.xaxis.set_ticks(ticks)
    ax.xaxis.set_ticks_position("bottom")
    ax.xaxis.set_tick_params(which='major', length=5 * scale)
    ax.xaxis.set_ticklabels(ticks,
                            fontproperties=textprops)
    ax.set_xlabel("TIME",  fontproperties=textprops)
    ax.set_frame_on(False)
    return ax
