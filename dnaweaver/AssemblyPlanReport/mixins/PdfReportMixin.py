import os
import datetime
import weasyprint
import pandas
from jinja2 import Environment, FileSystemLoader
from ..config import SETTINGS
from .PlotsMixin.give_quotes_html_locations import give_quotes_html_locations
from ...version import __version__ as version


class PdfReportMixin:
    def make_html_report(self):
        """Return a HTML version of the assembly report (later converted to PDF)."""

        def _anchor_span(op):
            return "<span id='anchor_%s'>%s</span>" % (op.id, op.id)

        operations = self.to_steps_list()

        def sorting_key(op):
            location = op.get("final_location", None)
            if location is None:
                location = (-1, 1)
            return (location[0], -location[1])

        operations = sorted(operations, key=sorting_key)

        give_quotes_html_locations(operations, len_sequence=len(self.plan.sequence))
        env = Environment(loader=FileSystemLoader(SETTINGS["template_path"]))
        template = env.get_template("report_template.html")
        orders_dataframe = pandas.DataFrame.from_records(
            [
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
                if self.sources[op.source].operation_type == "order"
            ],
            columns=[
                "ID",
                "Company",
                "Sequence",
                "Length",
                "Sequence",
                "Price",
                "Lead time",
                "Location",
            ],
        )

        pcrs_dataframe = pandas.DataFrame.from_records(
            [
                {
                    "ID": _anchor_span(op),
                    "Primers": ", ".join(
                        [_anchor_span(child) for child in op.get("assembly_plan", [])]
                    ),
                    "Source": op.source,
                    "Sequence": op.sequence,
                    "Length": len(op.sequence),
                    "Price": op.price,
                    "Lead time": op.lead_time,
                    "Deadline": op.get("deadline", None),
                    "Location": op.html_location,
                }
                for op in operations
                if self.sources[op.source].operation_type == "PCR"
            ],
            columns=[
                "ID",
                "Source",
                "infos",
                "Primers",
                "Sequence",
                "Length",
                "Price",
                "Lead time",
                "Location",
            ],
        )

        parts_dataframe = pandas.DataFrame.from_records(
            [
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
                if self.sources[op.source].operation_type == "library"
            ],
            columns=[
                "ID",
                "Source",
                "Sequence",
                "Length",
                "Price",
                "Lead time",
                "Location",
            ],
        )

        assembly_operations = [
            op
            for op in operations
            if self.sources[op.source].operation_type == "assembly"
        ]

        asm_dataframe = pandas.DataFrame.from_records(
            [
                {
                    "ID": _anchor_span(op),
                    "Station": op.source,
                    "Sequence": op.sequence,
                    "Length": len(op.sequence),
                    "Price": op.price,
                    "Lead time": op.lead_time,
                    "Deadline": op.get("deadline", None),
                    "Location": op.html_location,
                    "Primers": None,
                }
                for op in assembly_operations
            ],
            columns=[
                "ID",
                "Station",
                "infos",
                "Sequence",
                "Length",
                "Price",
                "Lead time",
                "Location",
            ],
        )

        cost_orders = orders_dataframe["Price"].sum()

        render_parameters = {
            "sequence": self.plan.sequence,
            "sequence_name": self.plan.get("sequence_name", "unknown"),
            "sequence_hash": "%02x" % hash(self.plan.sequence),
            "quote": self.plan,
            "date_submitted": datetime.datetime.now(),
            "n_orders": len(orders_dataframe),
            "cost_orders": cost_orders,
            "total_cost": self.plan.price,
            "orders": orders_dataframe,
            "asm_dataframe": asm_dataframe,
            "assembly_operations": assembly_operations,
            "user_name": self.plan.get("user_name", "unknown"),
            "version": version,
        }
        if len(pcrs_dataframe) > 0:
            render_parameters["pcr_reuses"] = pcrs_dataframe

        if len(parts_dataframe) > 0:
            render_parameters["parts_reuses"] = parts_dataframe

        html_out = template.render(render_parameters)

        return html_out

    def write_pdf_report(self, target):
        """Return a PDF version of the report with general infos and details of
        each intermediary constructs."""
        base_path = SETTINGS["template_path"]
        css_path = os.path.join(base_path, "static", "css")
        stylesheets = [
            os.path.join(css_path, "dnaweaver_report.css"),
            # os.path.join(css_path, "font-awesome.min.css")
        ]
        html_out = self.make_html_report()
        weasy_html = weasyprint.HTML(string=html_out, base_url=base_path)
        weasy_html.write_pdf(target, stylesheets=stylesheets)
