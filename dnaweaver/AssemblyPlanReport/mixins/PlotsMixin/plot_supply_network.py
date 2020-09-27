from collections import defaultdict
import textwrap
from ...config import SETTINGS
from .plot_tree_graph import plot_tree_graph
import matplotlib.font_manager as fm


def plot_supply_network(
    quote=None,
    main_source=None,
    edges=None,
    levels=None,
    ax=None,
    textprops=None,
    margin=None,
    scale=1.0,
    interlevel_shift=0,
):
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
      Properties of the text. If None, a nice opensans font is used.

    Returns
    -------

    elements_positions, ax
      Dictionary of element positions, matplotlib ax.
    """
    if main_source is not None:
        edges, levels = main_source.compute_supply_graph()
    if edges is None:
        edges = []
        levels = defaultdict(lambda *a: [])
        for source_name, infos in quote.sources.items():
            levels[infos._depth].append(source_name)
            for provider_name in infos.providers:
                edges.append((provider_name, source_name))
        levels = [levels[i] for i in range(max(levels) + 1)][::-1]
    else:
        sources_dict = {source.name: source for level in levels for source in level}

    fontawesome = fm.FontProperties(
        fname=SETTINGS["fontawesome-ttf-path"], size=18 * scale, family="sans-serif",
    )
    if textprops is None:
        textprops = fm.FontProperties(
            fname=SETTINGS["OpenSans-ttf-path"], size=12 * scale, family="sans-serif",
        )

    def draw_node(x, y, source_name, ax):
        if quote:
            source = quote.sources[source_name]
            symbol = source._report_fa_symbol
        else:
            source = sources_dict[source_name.name]
            symbol = source.report_fa_symbol

        ax.text(
            x,
            y,
            symbol,
            horizontalalignment="center",
            verticalalignment="center",
            fontproperties=fontawesome,
        )

        text = "\n".join(textwrap.wrap(source.name, 13))

        ax.text(
            x,
            y + 0.035 * scale ** 0.3,
            text,
            horizontalalignment="center",
            verticalalignment="bottom",
            fontproperties=textprops,
        )

    return plot_tree_graph(
        levels,
        edges,
        draw_node,
        ax=ax,
        scale=scale,
        interlevel_shift=interlevel_shift,
        margin=margin,
    )
