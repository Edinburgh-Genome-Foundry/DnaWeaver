"""Plotting functions for the results of DNA Weaver."""

from io import StringIO, BytesIO
from base64 import b64encode


def matplotlib_figure_to_file_string(fig, format="svg", **kwargs):
    """Return a string of the figure in the requested format."""

    if format == "pdf":
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
