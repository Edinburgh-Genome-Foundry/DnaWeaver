"""Plotting functions for the results of DNA Weaver"""

import textwrap
from copy import deepcopy
import itertools
from collections import defaultdict

import colorsys

from Bio import SeqIO

try:
    from StringIO import StringIO
    USE_BYTES = False
except ImportError:
    from io import StringIO, BytesIO
    USE_BYTES = True

from base64 import b64encode

try:
    from bokeh.plotting import figure, ColumnDataSource
    from bokeh.models import (Range1d, HoverTool, Legend)
    BOKEH_AVAILABLE = True
except ImportError:
    BOKEH_AVAILABLE = False
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cl
import matplotlib.cm as cm
import matplotlib.font_manager as fm
from matplotlib.path import Path
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

from dna_features_viewer import GraphicRecord

from .config import SETTINGS


def plot_assembly_blocks(quote, parts_offset=0, backend="matplotlib",
                         plot_top_assembly=True,
                         ax=None, edge_widths=None, legend=False):
    """Return a Matplotlib or Bokeh plot of the assembly tree of blocks.

    Parameters
    ----------



    parts_offset
      Offset applied so that consecutive blocks are not exactly on the same
      level. Can go from 0 (flat line of blocks) to e.g. 1

    backend
      One of "matplotlib" or "bokeh". Determines the library used to produce
      the graph. matplotlib can output PNG/SVG/etc. while Bokeh outputs an
      interactive HTML/JS graph

    plot_top_assembly
      Whether the top assembly (which is just one big block) should be plotted
      or not

    ax
      A Matplotlib Axes object. If no ax is provided, a new figure and ax are
      generated.

    edge_widths
      A dict {sourcename : width} indicating the widths of the rectangles based
      on the source name. Rectangles with a very small width will be edgeless.

    legend
      Whether the legend is included in the ax.

    """
    rectangles = []
    if edge_widths is None:
        edge_widths = {}
    if not hasattr(quote, 'tree'):
        print('lol')
        quote.compute_full_assembly_tree()
    tree = deepcopy(quote.tree)

    def rec(_quote, depth=0, xstart=0, xend=None):
        children = _quote.pop("assembly_plan")
        #print _quote.final_location
        if _quote.final_location is not None:
            left, right = _quote.final_location
        else:
            left, right = xstart, xend

        if children is None:
            children = []
        if plot_top_assembly or (depth > 0):
            source = quote.sources[_quote.source]
            rectangles.append({
                "top": -depth + .5 + (not plot_top_assembly),
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
                "color": source._report_color
            })
        for child in children:
            rec(child, depth=depth + 1,
                xstart=xstart + child.segment_start,
                xend=xstart + child.segment_end)
    rec(tree, xend=len(tree.sequence))
    tops = list(set([r["top"] for r in rectangles]))
    offsets = {top: 0 for top in tops}

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
        if parts_offset:
            offset = offsets[rectangle["top"]]
            offsets[rectangle["top"]] = (offsets[rectangle["top"]] + 1) % 2
            rectangle["top"] += parts_offset * offset
            rectangle["bottom"] += parts_offset * offset

    if backend == "bokeh":
        return _bokeh_plot_assembly_blocks(rectangles, guides, tops)
    elif backend == "matplotlib":
        if hasattr(quote, "genbank") and (quote.genbank is not None):
            record = SeqIO.read(quote.genbank, "genbank")
        else:
            record = None
        return _matplotlib_plot_assembly_blocks(
            rectangles, guides, tops, ax, record=record, legend=legend)
    else:
        raise ValueError("Backend should be bokeh or matplotlib")


def _bokeh_plot_assembly_blocks(rectangles, guides, tops):
    """Plot the block rectangles in a bokeh plot"""
    hover = HoverTool(
        tooltips="""
        <div>
          <i class="fa fa-@fa_symbol"></i>
          <b>@operation_type from @source </b><br/>
          @length bp, @price $, @lead_time days
          </div>
    """,
        names=["parts"]
    )

    L = max([r["right"] for r in rectangles])
    p = figure(tools=[hover, "xwheel_zoom,xpan,reset"],
               responsive=True,
               plot_height=300,
               x_range=Range1d(int(-0.05 * L), int(1.05 * L)),
               y_range=Range1d(min(tops) - 1, 1),
               logo=None, toolbar_location="above",
               title="Assembly Plan")
    source_guide = ColumnDataSource(pd.DataFrame.from_records(guides))
    p.quad(name="guides",
           top="top", bottom="bottom", left="left",
           right="right", color="color", source=source_guide, alpha=0.5)
    dataframe = pd.DataFrame.from_records(rectangles)
    legend_items = []
    for sourcename, df in dataframe.groupby("source"):
        source = ColumnDataSource(df)
        quads = p.quad(name="parts", top="top", bottom="bottom", left="left",
                       right="right", color="color", source=source,
                       line_width=1, line_color="black", alpha=1)
        legend_items.append([sourcename, [quads]])

    p.yaxis.visible = False
    p.outline_line_color = None
    p.grid.grid_line_color = None
    p.add_layout(Legend(items=legend_items, location=(0, -30)), 'right')
    return p


def _matplotlib_plot_assembly_blocks(rectangles, guides, tops, ax=None,
                                     fig_height="auto", fig_width=8,
                                     legend=False, textprops=None,
                                     record=None):
    """Plot the assembly block rectangles using matplotlib."""

    if record is not None:
        gs = gridspec.GridSpec(4, 1)
        fig = plt.figure(figsize=(7,7))
        record_ax = fig.add_subplot(gs[0])
        blocks_ax = fig.add_subplot(gs[1:])
        grecord = GraphicRecord.from_biopython_record(
            record, fun_color=lambda a: "#aabbff"
        )
        grecord.plot(record_ax)
        res = plot_assembly_blocks(rectangles, guides, tops, ax=blocks_ax,
                                   fig_height="auto", fig_width=8,
                                   legend=False, textprops=None,
                                   record=None)
        record_ax.set_xlim(blocks_ax.get_xlim())
        return res

    L = max([r["right"] for r in rectangles])

    if textprops is None:
        textprops = fm.FontProperties(
            fname=SETTINGS["OpenSans-ttf-path"],
            size=12, family="sans-serif"
        )

    if ax is None:
        if fig_height == "auto":
            fig_height = 8
        fig, ax = plt.subplots(1, figsize=(8, 8))
    ax.set_xlim((-0.05 * L, 1.05 * L))
    ax.set_ylim((min(tops) - 1, 1.5))

    for g in guides:
        patch = mpatches.Rectangle((g["left"], g["bottom"]),
                                   g["right"] - g["left"],
                                   g["top"] - g["bottom"],
                                   color=g["color"],
                                   ec="none")
        ax.add_patch(patch)

    seen_sources = set([])
    legend_handles = []
    sorted_rectangles = sorted(rectangles,
                               key=lambda r: (-r["bottom"], r["left"]))
    for g in sorted_rectangles:
        if g["source"] not in seen_sources:
            legend_patch = mpatches.Patch(facecolor=g["color"],
                                          label=g["source"],
                                          edgecolor="k",
                                          linewidth=1.0)
            legend_handles.append(legend_patch)
            seen_sources.add(g["source"])
            plt.legend(handles=[legend_patch])
        # linewidth = g["line_width"]
        width = g["right"] - g["left"]
        line_width = 1.0 if (1.0 * width / L) > .002 else 0
        patch = mpatches.Rectangle((g["left"], g["bottom"]),
                                   width,
                                   g["top"] - g["bottom"],
                                   facecolor=g["color"],
                                   edgecolor="k",
                                   linewidth=line_width)
        ax.add_patch(patch)
    ax.axis("off")

    if legend:
        if ax.legend_ is not None:
            ax.legend_.remove()
        legend = ax.legend(handles=legend_handles,
                           frameon=False, ncol=2, loc=2,
                           bbox_to_anchor=(0.0, -.05))
        ltext = legend.get_texts()
        plt.setp(ltext, fontproperties=textprops)
        ax.figure.subplots_adjust(hspace=.0)

    return ax, legend_handles


def plot_assembly_timeline(quote, deadline=None, ax=None,
                           rectangle_color="#bbbbff", scale=1.0,
                           backend="matplotlib"):
    """Make a Gantt-like chart of the assemblies."""
    if backend == "matplotlib":
        return _matplotlib_plot_assembly_timeline(
            quote, deadline=deadline, ax=ax, rectangle_color=rectangle_color,
            scale=scale
        )


def _matplotlib_plot_assembly_timeline(quote, deadline=None, ax=None,
                                       rectangle_color="#bbbbff", scale=1.0):
    """Plot the assembly timeline with Matplotlib."""
    if deadline is None:
        deadline = quote.lead_time
    if not hasattr(quote, 'tree'):
        quote.compute_full_assembly_tree()
    quote.propagate_deadline(deadline)
    assemblies_list = quote.tree_as_list()
    assemblies_list = sorted(assemblies_list, key=lambda a: a.id)

    textprops = fm.FontProperties(
        fname=SETTINGS["OpenSans-ttf-path"],
        size=12 * scale, family="sans-serif"
    )
    fontawesome = fm.FontProperties(
        fname=SETTINGS["fontawesome-ttf-path"],
        size=12 * scale, family="sans-serif"
    )
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(16 * scale,
                                           len(assemblies_list) * 0.3 * scale))
    positions = {}

    if textprops is None:
        textprops = fm.FontProperties(
            fname=SETTINGS["OpenSans-ttf-path"],
            size=12, family="sans-serif"
        )

    for i, assembly in enumerate(assemblies_list):
        source = assembly.source
        color = "k"
        if source.__class__.__name__ == "CommercialDnaOffer":
            color = "#888888"
        x2 = assembly.deadline
        x1 = assembly.deadline - assembly.step_duration
        positions[assembly.id] = ([x1, i], [x2, i])
        patch = mpatches.Rectangle(
            (assembly.deadline - assembly.step_duration, i - .2),
            assembly.step_duration, .4, color=rectangle_color,
            ec="none"
        )

        ax.plot([0, assembly.deadline], [i, i], c="k", alpha=0.3,
                lw=.5 * scale, zorder=-1000, ls="-")
        ax.add_patch(patch)

        ax.text(-1.5, i, assembly.id, horizontalalignment="right",
                verticalalignment="center", fontproperties=textprops,
                color=color)
        ax.text(-.75, i, source.report_fa_symbol, horizontalalignment="center",
                verticalalignment="center", fontproperties=fontawesome,
                color=color)

    for assembly in assemblies_list:
        xa, ya = positions[assembly.id][0]
        if assembly.assembly_plan is not None:
            for segment, child in assembly.assembly_plan.items():
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
    ax.set_xlabel("TIME", fontproperties=textprops)
    ax.set_frame_on(False)
    return ax


def plot_tree_graph(levels, edges, draw_node, elements_positions=None,
                    ax=None, width_factor=2.5, height_factor=2, scale=1.0,
                    edge_left_space=0.015, edge_right_space=0.015,
                    interlevel_shift=0, margin=None, **txt_kw):
    """General function for plotting tree graphs.

    Parameters
    ----------

    levels
      A list of lists of nodes grouped by "level", i.e distance to the in the
      graph to the level 0. levels will be displayed on a same column.

    edges
      List of nodes pairs (source node, target node).

    draw_node
      A function f(x , y , node, ax, **kw) which draws something related to the
      node at the position x,y on the given Matplotlib ax.

    ax
      The matplotlib ax to use. If none is provided, a new ax is generated.

    Returns
    -------
    elements_positions, ax
      Dictionary of elements positions, matplotlib ax.

    Examples:
    ---------

    >>> def draw_node(x,y, node, ax):
        ax.text(x,y, node)
    >>> plot_tree_graph(levels=[["A","B","C"], ["D,E"], ["F"]],
                        edges=[("A","D"),("B","D"),("C","E")
                               ("D","F"),("E","F")],
                        draw_node = draw_node,)



    """
    levels_dict = {
        element: level
        for level, elements in enumerate(levels)
        for element in elements
    }
    if elements_positions is None:
        elements_positions = {}
        for lvl, elements in enumerate(levels):
            yy = np.linspace(0, 1, len(elements) + 2)[1:-1]
            yy += interlevel_shift * (1-2*(lvl % 2))
            x = 1.0 * (1 + lvl) / (len(levels) + 1)
            for y, element in zip(yy, elements):
                elements_positions[element] = (x, y)

    if ax is None:
        width = width_factor * len(levels) * scale
        height = height_factor * max([len(lvl) for lvl in levels]) * scale
        fig, ax = plt.subplots(1, figsize=(width, height))

    for element, (x, y) in elements_positions.items():
        draw_node(x, y, element, ax, **txt_kw)

    y_spans = [
        elements_positions[elements[1]][1] - elements_positions[elements[0]][1]
        for elements in levels
        if len(elements) > 1
    ]

    delta_y = 0.5*min(y_spans) if y_spans != [] else 0

    for el1, el2 in edges:
        x1, y1 = elements_positions[el1]
        x2, y2 = elements_positions[el2]
        x1 += edge_left_space * np.sqrt(scale)
        x2 += -edge_right_space * np.sqrt(scale)
        if ((levels_dict[el2] - levels_dict[el1]) > 1) and (y1 == y2):
            patch = mpatches.PathPatch(
                Path([(x1, y1), (0.5 * x2 + 0.5 * x1, y1-delta_y),
                      (0.5 * x2 + 0.5 * x1, y2-delta_y), (x2, y2)],
                     [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                facecolor='none', lw=1 * scale,
                zorder=-1000
            )

        else:
            patch = mpatches.PathPatch(
                Path([(x1, y1), (0.9 * x2 + 0.1 * x1, y1),
                      (0.1 * x2 + 0.9 * x1, y2), (x2, y2)],
                     [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                facecolor='none', lw=1 * scale,
                zorder=-1000
            )

        ax.add_patch(patch)

    ax.axis("off")
    if margin is not None:
        xx, yy = [np.array(e) for e in zip(*elements_positions.values())]
        xmin, xmax = xx.min(), xx.max()
        dx = margin * (xmax - xmin)
        ymin, ymax = yy.min(), yy.max()
        dy = margin * (ymax - ymin)
        d = max(dx, dy)
        ax.set_xlim(xmin - d, xmax + d)
        ax.set_ylim(ymin - d, ymax + d)

    return elements_positions, ax


def plot_supply_graph(quote, ax=None, textprops=None, margin=None,
                      scale=1.0, interlevel_shift=0):
    """Plot the supply graph of all sources related to the provided source.

    Examples:
    ---------

    >>> quote.sources = chunks_assembly_station.compute_dict_supply_graph()

    Parameters
    ----------

    quote.sources
      A dictionary {"source_1_name": {"depth": 1, "infos": {...}},
                    "source_2_name": {"depth": 2, "infos": {...}},
                    etc.}
    ax
      A matplotlib ax. If None provided, one is created (and returned at the
      end)

    textprops
      Properties of the text. If none, a nice opensans font is used.
    
    Returns
    -------
    elements_positions, ax
      Dictionary of elements positions, matplotlib ax.


    """
    edges = []
    levels = defaultdict(lambda *a: [])
    for source_name, infos in quote.sources.items():
        levels[infos._depth].append(source_name)
        for provider_name in infos.providers:
            edges.append((provider_name, source_name))
    levels = [levels[i] for i in range(max(levels) + 1)][::-1]

    fontawesome = fm.FontProperties(
        fname=SETTINGS["fontawesome-ttf-path"],
        size=18 * scale, family="sans-serif"
    )
    if textprops is None:
        textprops = fm.FontProperties(
            fname=SETTINGS["OpenSans-ttf-path"],
            size=12 * scale, family="sans-serif"
        )

    def draw_node(x, y, source_name, ax):
        source_infos = quote.sources[source_name]
        symbol = source_infos._report_fa_symbol
        ax.text(x, y, symbol, horizontalalignment="center",
                verticalalignment="center", fontproperties=fontawesome)

        text = "\n".join(textwrap.wrap(source_infos.name, 13))

        ax.text(x, y + 0.035 * scale**.3, text,
                horizontalalignment="center",
                verticalalignment="bottom", fontproperties=textprops,
                )
    return plot_tree_graph(levels, edges, draw_node, ax=ax, scale=scale,
                           interlevel_shift=interlevel_shift, margin=margin)


def plot_assembly_graph(quote, ax=None, margin=None, textprops=None,
                        scale=1.0):
    """Plot the complete assembly graph
    
    Returns
    -------
    elements_positions, ax
      Dictionary of elements positions, matplotlib ax.
    """

    nodes_dict = {}
    levels = defaultdict(lambda *a: [])
    edges = []
    if not hasattr(quote, 'tree'):
        quote.compute_full_assembly_tree()
    tree = deepcopy(quote.tree)

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
        size=13 * scale, family="sans-serif"
    )
    if textprops is None:
        textprops = fm.FontProperties(
            fname=SETTINGS["OpenSans-ttf-path"],
            size=12 * scale, family="sans-serif"
        )

    def draw_node(x, y, node_id, ax):
        node = nodes_dict[node_id]
        icon = quote.sources[node.source]._report_fa_symbol
        ax.text(x, y, node_id, horizontalalignment="left",
                verticalalignment="center",
                fontproperties=textprops)
        ax.text(x - 0.01 * np.sqrt(scale), y, icon,
                horizontalalignment="right",
                verticalalignment="center",
                fontproperties=fontawesome)

    all_elements = sorted(sum(levels, []))
    ypos = {
        el: 1.0 * (i + 1) / (len(all_elements) + 2)
        for i, el in enumerate(all_elements)
    }
    for el in all_elements:
        children = [
            e2
            for (e2, e1) in edges
            if e1 == el
        ]
        if children != []:
            ypos[el] = 1.0 * sum(ypos[e] for e in children) / len(children)

    xpos = {
        el: 1.0 * (1 + x) / (len(levels) + 1)
        for x, elements in enumerate(levels)
        for el in elements
    }

    elements_positions = {el: (xpos[el], ypos[el]) for el in all_elements}
    return plot_tree_graph(levels, edges, draw_node,
                           elements_positions=elements_positions, ax=ax,
                           edge_left_space=0.06,
                           edge_right_space=0.03,
                           margin=margin,
                           height_factor=0.40, width_factor=5.5, scale=scale)


def matplotlib_figure_to_file_string(fig, format="svg", **kwargs):
    """Return a string of the figure in the requested format."""

    if (format == "pdf") and USE_BYTES:
        output = BytesIO()
    else:
        output = StringIO()

    fig.savefig(output, format=format, **kwargs)
    return output.getvalue()


def matplotlib_figure_to_svg_base64_data(fig, **kwargs):
    """Return a string of the form 'data:image/svg+xml;base64,XXX' where XXX
    is the base64-encoded svg version of the figure."""
    svg_txt = matplotlib_figure_to_file_string(fig, format="svg", **kwargs)
    svg_txt = "\n".join(svg_txt.split("\n")[4:])
    svg_txt = "".join(svg_txt.split("\n"))
    try:
        return "data:image/svg+xml;base64," + b64encode(svg_txt)
    except:
        content = b64encode(svg_txt.encode("ascii"))
        result = (b"data:image/svg+xml;base64," + content).decode("utf-8")
        return str(result)


def give_quotes_html_locations(quotes, len_sequence, ax=None):
    """Quickly adds an html-location field containing an HTML image of the
    location of each sub-quote in the final sequence.

    The image is a line representing the full sequence and a rectangle
    indicating the position of the subquote in the the full sequence.

    The advantage of processing all quotes at once is that the
    figure is drawn once, then only the position of the rectangle and the
    title are modified.
    """
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(1.5, 0.2), facecolor=None)
    ax.set_frame_on(False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.plot([0, len_sequence], [0, 0], c="k")
    rect = mpatches.Rectangle((0, -0.5), 0, 1, color="#ddddff",
                              ec="k", alpha=1, zorder=1000,
                              transform=ax.transData)
    ax.add_patch(rect)
    ax.set_ylim((-1, 1))
    ax.set_xlim(0, len_sequence)
    text = ax.text(0.5, 0.65, "none", horizontalalignment="center", fontsize=4)
    for q in quotes:
        location = q.final_location
        if location is None:
            q.html_location = ''
        else:
            rect.set_bounds(location[0], -.5, location[1] - location[0], 1.0)
            text.set_text("%d-%d" % (location[0], location[1]))
            svg_data = matplotlib_figure_to_svg_base64_data(ax.figure,
                                                            transparent=True)
            q.html_location = "<img src='%s' />" % svg_data
    plt.close(ax.figure)


def hls_to_hex(hue, luminance, saturation):
    """Return (R,G,B) equivalent of a hue/staturation/value color."""
    return cl.rgb2hex(colorsys.hls_to_rgb(hue, luminance, saturation))

def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (int(red), int(green), int(blue))

def autocolor_quote_sources(quote, hues=(0.635, 0.047, 0.117),
                            saturations=(0.9, 0.7, 0.5, 0.3),
                            min_lum=0.2, max_lum=0.8):
    """Auto-add a `_report_color` field to the sources in in quote.sources.

    Sources at the same depth share the same luminance
    """

    colors = itertools.cycle([
        rgb_to_hex(*[
            255*e**0.4
            for e in cm.Paired(0.13 * i % 1.0)
        ][:3])
        for i in range(30)
    ])
    for name, source in sorted(quote.sources.items()):
        color = next(colors)
        source._report_color = color

def ax_to_pdf(ax):
    return lambda fh: ax.figure.savefig(fh, format="pdf", bbox_inches="tight")

def plot_decomposition_graph(graph, nodes_color="#6886b7", weight='weight',
                             colormap='jet', ax=None, edge_width=1,
                             edges_step=1, edge_alpha=0.2, figsize=(8, 8)):
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
    normalized_weights = (255 * weights / weights.max()).astype('uint8')
    colormap = cm.__dict__[colormap]
    colors = colormap(normalized_weights, alpha=edge_alpha)
    if ax is None:
        fig, ax = plt.subplots(1, figsize=figsize)
    ax.axis("off")
#     print (list(colors))
    for (start, end, w), color in zip(edges, colors):
        xc = 0.5 * (start + end)
        half = 0.5 * abs(end - start)
        ax.add_patch(mpatches.Arc((xc, 0), 2 * half, 2 * half,
                                  theta1=0 ,theta2=180, facecolor='none',
                                  ls='-', edgecolor=color, linewidth=1))
    ax.plot(xx, [0 for x in xx], marker="o", c=nodes_color)
    ax.set_aspect("equal")
    ax.set_xlim(-1, L + 1)
    ax.set_ylim(-1, max_segment_length / 2 + 1)

def plot_cost_profiles(cost_profiles):
    """Plot cost profiles generated in OptimizeManufacturability"""
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    L = len(cost_profiles['count'])
    diffs = cost_profiles['diffs']
    r = (L - len(diffs)) / 2
    xx = range(L)
    ax1.plot(r + np.arange(len(diffs)), diffs, c='orange')
    ymax = diffs.max()
    if ymax == 0:
        ymax = 1
    ax1.set_ylim(bottom=-0.005, top = 1.5 * ymax)
    ax1b = ax1.twinx()
    ax1b.plot(xx, cost_profiles['count'], c='k', alpha=0.3)
    ax1b.set_ylim(bottom=0, top = 1.5 * cost_profiles['count'].max())
    ax1.set_xlim(0, L)
    
    ax2.plot(xx, cost_profiles['mean'], c='b', lw=2)
    ax2.fill_between(
        xx,
        cost_profiles['min'],
        cost_profiles['max'],
        facecolor='b',
        alpha=0.2
    )
    ax2.set_ylim(bottom=0, top = 1.5 * cost_profiles['max'].max())
    return ax1, ax2