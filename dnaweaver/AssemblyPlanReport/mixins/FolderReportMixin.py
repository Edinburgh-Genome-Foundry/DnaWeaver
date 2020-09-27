import os
import flametree
import pandas
from ..config import SETTINGS


class FolderReportMixin:
    def write_full_report(self, target="@memory"):
        """Generate an extensive, multifile report as either a folder or a zip."""

        def write_ax_as_pdf(ax, target):
            ax.figure.savefig(target, format="pdf", bbox_inches="tight")

        root = flametree.file_tree(target, replace=True)
        self.write_pdf_report(root._file("assembly_report.pdf").open("wb"))
        self.write_sequences_spreadsheet(root._file("sequences.csv").open("w"))
        root._file("README.txt").write(_get_folder_readme_content())
        self.write_figures_folder(root._dir("figures"))
        self.write_all_sequence_records(root._dir("genbank"))
        return root._close()

    def write_sequences_spreadsheet(quote, target):

        quotes_list = quote.to_steps_list()
        sequences_dataframe = pandas.DataFrame.from_records(
            [
                {
                    "ID": _quote.id,
                    "Type": _quote.get("operation_type", None),
                    "Source": _quote.source,
                    "SequenceLength": len(_quote.sequence),
                    "Sequence": _quote.sequence,
                }
                for _quote in quotes_list
            ],
            columns=["ID", "Type", "Source", "SequenceLength", "Sequence"],
        ).sort_values(by=["Type", "Source", "ID"])

        sequences_dataframe.to_csv(target, sep=";", index=False)


def _get_folder_readme_content():
    """General function to get the content of the FOLDER_README.txt template."""
    readme_path = os.path.join(SETTINGS["template_path"], "FOLDER_README.txt")
    with open(readme_path, "r") as f:
        readme_content = f.read()
    return readme_content
