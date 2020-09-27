import matplotlib.pyplot as plt
from .AssemblyBlocksMixin import AssemblyBlocksMixin
from .AssemblyGraphMixin import AssemblyGraphMixin
from .ColorsMixin import ColorsMixin
from .TimelineMixin import TimelineMixin

from .plot_supply_network import plot_supply_network

import warnings

warnings.filterwarnings("ignore", message="Glyph 112 missing from current font.")
warnings.filterwarnings("ignore", message="Glyph 108 missing from current font.")


class PlotsMixin(AssemblyBlocksMixin, AssemblyGraphMixin, ColorsMixin, TimelineMixin):
    def write_figures_folder(self, figures_folder):
        def write_ax_as_pdf(ax, target):
            ax.figure.savefig(target.open("wb"), format="pdf", bbox_inches="tight")
            plt.close(ax.figure)

        pos, ax = self.plot_supply_network()
        write_ax_as_pdf(ax, figures_folder._file("supply_network.pdf"))
        plt.close(ax.figure)

        pos, ax = self.plot_assembly_graph(ax=None, textprops=None)
        write_ax_as_pdf(ax, figures_folder._file("assembly_graph.pdf"))

        _, ax_assembly_blocks = plt.subplots(1, figsize=(7, 4))
        assembly_blocks_ax, lg = self.plot_assembly_blocks(
            plot_top_assembly=False,
            ax=ax_assembly_blocks,
            parts_offset=0.1,
            legend=True,
        )
        assembly_blocks_ax.figure.subplots_adjust(bottom=0.3)
        write_ax_as_pdf(assembly_blocks_ax, figures_folder._file("assembly_blocks.pdf"))

    def plot_supply_network(self, ax=None):
        """Plot the supply network (see plot_supply_network.plot_supply_network
        for more options).

        Returns
        -------

        elements_positions, ax
          Dictionary of elements positions, matplotlib ax.
        """

        return plot_supply_network(self, ax=ax)
