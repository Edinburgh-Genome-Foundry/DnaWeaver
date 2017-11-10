
import weasyprint
import matplotlib.pyplot as plt
from copy import deepcopy
import datetime
import os
import flametree

from jinja2 import Environment, FileSystemLoader
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

try:  # Python2 vs. Python3 ways of dealing with file-like strings
    from StringIO import StringIO
    USE_BYTES = False
except ImportError:
    from io import StringIO, BytesIO
    USE_BYTES = True

import pandas as pd

from .config import SETTINGS
from .plotting import (plot_assembly_blocks,
                       plot_supply_graph,
                       plot_assembly_graph,
                       give_quotes_html_locations,
                       ax_to_pdf)


pd.set_option('display.max_colwidth', -1)


def quote_to_SeqRecord(quote, record=None, record_id=None):
    """Return a Biopython seqrecord of the quote.

    >>> record = to_SeqRecord(solution)
    >>> # Let's plot with DnaVu:
    >>> from dnavu import create_record_plot
    >>> from bokeh.io import output_file, show
    >>> output_file("view.html")
    >>> plot = create_record_plot(record)
    >>> show(plot)
    """
    if record_id is None:
        record_id = quote.id
    if record is None:
        record = SeqRecord(Seq(quote.sequence, DNAAlphabet()), id=record_id)
    else:
        record = deepcopy(record)

    if quote.assembly_plan is not None:
        features = [
            SeqFeature(
                FeatureLocation(q.segment_start, q.segment_end, 1),
                type="Feature",
                qualifiers={
                    "name": q.id,
                    "source": q.source,
                    "price": q.price,
                    "lead_time": q.lead_time
                }
            )
            for q in quote.assembly_plan
        ]
        record.features = features + record.features
    return record


def quote_to_genbank(quote, filename=None, filehandle=None,
                     record=None, record_id=None):
    record = quote_to_SeqRecord(quote, record=record, record_id=record_id)
    if filename is not None:
        with open(filename, "w+") as f:
            SeqIO.write(record, f, "genbank")
    else:
        output = StringIO()
        SeqIO.write(record, output, "genbank")
        return output.getvalue()


def make_html_report(quote):
    """Return a HTML version of the assembly report (later converted to PDF)"""

    def _anchor_span(op):
        return "<span id='anchor_%s'>%s</span>" % (op.id, op.id)

    operations = quote.to_quotes_list()

    def sorting_key(op):
        location = op.get("final_location", None)
        if location is None:
            location = (-1, 1)
        return (location[0], -location[1])

    operations = sorted(operations, key=sorting_key)

    give_quotes_html_locations(
        operations, len_sequence=len(quote.tree.sequence))

    path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        "templates")

    env = Environment(loader=FileSystemLoader(path))
    template = env.get_template("dnaweaver_jinja.html")

    orders_dataframe = pd.DataFrame.from_records([
        {
            "ID": _anchor_span(op),
            "Company": op.source,
            "Sequence": op.sequence,
            "Length": len(op.sequence),
            "Price": op.price,
            "Lead time": op.lead_time,
            "Deadline": op.get("deadline", None),
            "Location": op.html_location,
        }
        for op in operations
        if quote.sources[op.source].operation_type == "order"
    ])

    pcrs_dataframe = pd.DataFrame.from_records([
        {
            "ID": _anchor_span(op),
            "Primers": ", ".join([
                _anchor_span(child)
                for child in op.get("assembly_plan", [])
            ]),
            "Source": op.source,
            "Sequence": op.sequence,
            "Length": len(op.sequence),
            "Price": op.price,
            "Lead time": op.lead_time,
            "Deadline": op.get("deadline", None),
            "Location": op.html_location,
        }
        for op in operations
        if quote.sources[op.source].operation_type == "PCR"
    ],
        columns=["ID", "Source", "infos", "Primers",
                 "Sequence", "Length", "Price", "Lead time", "Location"]
    )

    parts_dataframe = pd.DataFrame.from_records([
        {
            "ID": _anchor_span(op),
            "Source": op.source,
            "Sequence": op.sequence,
            "Length": len(op.sequence),
            "Price": op.price,
            "Lead time": op.lead_time,
            "Deadline": op.get("deadline", None),
            "Location": op.html_location,
        }
        for op in operations
        if quote.sources[op.source].operation_type == "library"
    ],
        columns=["ID", "Source", "Sequence", "Length", "Price",
                 "Lead time", "Location"]
    )

    assembly_operations = [
        op
        for op in operations
        if quote.sources[op.source].operation_type == "assembly"
    ]

    asm_dataframe = pd.DataFrame.from_records([
        {
            "ID": _anchor_span(op),
            "Station": op.source,
            "Sequence": op.sequence,
            "Length": len(op.sequence),
            "Price": op.price,
            "Lead time": op.lead_time,
            "Deadline": op.get("deadline", None),
            "Location": op.html_location,
            "Primers": None
        }
        for op in assembly_operations
    ],
        columns=["ID", "Station", "infos", "Sequence",
                 "Length", "Price", "Lead time", "Location"]
    )

    cost_orders = orders_dataframe["Price"].sum()

    render_parameters = {
        "sequence": quote.tree.sequence,
        "sequence_name": quote.tree.get("sequence_name", "unknown"),
        "sequence_hash": "%02x" % hash(quote.tree.sequence),
        "quote": quote.tree,
        "date_submitted": datetime.datetime.now(),
        "n_orders": len(orders_dataframe),
        "cost_orders": cost_orders,
        "total_cost": quote.tree.price,
        "orders": orders_dataframe,
        "asm_dataframe": asm_dataframe,
        "assembly_operations": assembly_operations,
        "user_name": quote.tree.get("user_name", "unknown")
    }
    if len(pcrs_dataframe) > 0:
        render_parameters["pcr_reuses"] = pcrs_dataframe

    if len(parts_dataframe) > 0:
        render_parameters["parts_reuses"] = parts_dataframe

    html_out = template.render(render_parameters)

    return html_out


def make_pdf_report(quote, filename=None):
    """Return a PDF version of the report with general infos and details of
    each intermediary constructs."""
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)))
    base_path = os.path.join(path, "templates")
    css_path = os.path.join(base_path, "static", "css")
    stylesheets = [
        os.path.join(css_path, "dnaweaver_report.css"),
        # os.path.join(css_path, "font-awesome.min.css")
    ]
    html_out = make_html_report(quote)
    weasy_html = weasyprint.HTML(string=html_out, base_url=base_path)
    if filename is not None:
        write_mode = "wb" if USE_BYTES else "w+"
        with open(filename, write_mode) as f:
            weasy_html.write_pdf(f, stylesheets=stylesheets)
    else:
        output = BytesIO() if USE_BYTES else StringIO()
        weasy_html.write_pdf(output, stylesheets=stylesheets)
        return output.getvalue()


def _get_folder_readme_content():
    """General function to get the content of the FOLDER_README.txt template"""
    readme_path = os.path.join(SETTINGS["template_path"], "FOLDER_README.txt")
    with open(readme_path, "r") as f:
        readme_content = f.read()
    return readme_content


def make_spreadsheet_sequences_report(quote, filename=None):

    quotes_list = quote.to_quotes_list()
    sequences_dataframe = pd.DataFrame.from_records(
        [
            {
                "ID": _quote.id,
                "Type": _quote.get("operation_type", None),
                "Source": _quote.source,
                "SequenceLength": len(_quote.sequence),
                "Sequence": _quote.sequence
            }
            for _quote in quotes_list
        ],
        columns=["ID", "Type", "Source", "SequenceLength", "Sequence"]
    ).sort_values(by=["Type", "Source", "ID"])

    result = sequences_dataframe.to_csv(sep=";", index=False)
    if filename is not None:
        with open(filename, "w+") as f:
            f.write(result)
    else:
        return result


def make_folder_report(quote, target="@memory"):
    """Generate an extensive, multifile report as either a folder or a zip."""

    def write_ax_as_pdf(ax, target):
        ax.figure.savefig(target, format="pdf", bbox_inches="tight")

    with flametree.file_tree(target) as root:
        root._file("assembly_report.pdf").write(
            make_pdf_report(quote))

        root._file('sequences.csv').write(
            make_spreadsheet_sequences_report(quote))

        root._file("README.txt").write(
            _get_folder_readme_content())

        # FIGURES

        figures = root._dir("figures")

        ax = plot_supply_graph(quote)
        write_ax_as_pdf(ax, figures._file("supply_network.pdf").open("wb"))
        plt.close(ax.figure)

        ax = plot_assembly_graph(quote, ax=None, textprops=None)
        write_ax_as_pdf(ax, figures._file("assembly_graph.pdf").open("wb"))
        plt.close(ax.figure)

        fig, ax_assembly_blocks = plt.subplots(1, figsize=(7, 4))
        assembly_blocks_ax, lg = plot_assembly_blocks(
            quote, backend="matplotlib", plot_top_assembly=False,
            ax=ax_assembly_blocks, parts_offset=0.1, legend=True)
        assembly_blocks_ax.figure.subplots_adjust(bottom=0.3)
        write_ax_as_pdf(assembly_blocks_ax,
                        figures._file("assembly_blocks.pdf").open("wb"))
        plt.close(fig)

        # GENBANKS

        genbank = root._dir("genbank")
        for q in quote.to_quotes_list():
            genbank._file(q.id + ".gb").write(quote_to_genbank(q))

    return root._close()
