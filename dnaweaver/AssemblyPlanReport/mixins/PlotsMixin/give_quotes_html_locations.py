import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from .matplotlib_export import matplotlib_figure_to_svg_base64_data


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
        _fig, ax = plt.subplots(1, figsize=(1.5, 0.2), facecolor=None)
    ax.set_frame_on(False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.plot([0, len_sequence], [0, 0], c="k")
    rect = mpatches.Rectangle(
        (0, -0.5),
        0,
        1,
        color="#ddddff",
        ec="k",
        alpha=1,
        zorder=1000,
        transform=ax.transData,
    )
    ax.add_patch(rect)
    ax.set_ylim((-1, 1))
    ax.set_xlim(0, len_sequence)
    text = ax.text(0.5, 0.65, "none", horizontalalignment="center", fontsize=4)
    for q in quotes:
        location = q.final_location
        if location is None:
            q.html_location = ""
        else:
            rect.set_bounds(location[0], -0.5, location[1] - location[0], 1.0)
            text.set_text("%d-%d" % (location[0], location[1]))
            svg_data = matplotlib_figure_to_svg_base64_data(ax.figure, transparent=True)
            q.html_location = "<img src='%s' />" % svg_data
    plt.close(ax.figure)
