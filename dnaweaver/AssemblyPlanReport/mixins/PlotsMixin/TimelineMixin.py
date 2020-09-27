import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

from ...config import SETTINGS


class TimelineMixin:
    def plot_assembly_timeline(
        self, deadline=None, ax=None, rectangle_color="#bbbbff", scale=1.0
    ):
        """Plot the assembly timeline with Matplotlib."""
        if deadline is None:
            deadline = self.plan.lead_time
        self.propagate_deadline(deadline)
        assemblies_list = self.tree_as_list()
        assemblies_list = sorted(assemblies_list, key=lambda a: a.id)

        textprops = fm.FontProperties(
            fname=SETTINGS["OpenSans-ttf-path"], size=12 * scale, family="sans-serif",
        )
        fontawesome = fm.FontProperties(
            fname=SETTINGS["fontawesome-ttf-path"],
            size=12 * scale,
            family="sans-serif",
        )
        if ax is None:
            _fig, ax = plt.subplots(
                1, figsize=(16 * scale, len(assemblies_list) * 0.3 * scale)
            )
        positions = {}

        if textprops is None:
            textprops = fm.FontProperties(
                fname=SETTINGS["OpenSans-ttf-path"], size=12, family="sans-serif"
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
                (assembly.deadline - assembly.step_duration, i - 0.2),
                assembly.step_duration,
                0.4,
                color=rectangle_color,
                ec="none",
            )

            ax.plot(
                [0, assembly.deadline],
                [i, i],
                c="k",
                alpha=0.3,
                lw=0.5 * scale,
                zorder=-1000,
                ls="-",
            )
            ax.add_patch(patch)

            ax.text(
                -1.5,
                i,
                assembly.id,
                horizontalalignment="right",
                verticalalignment="center",
                fontproperties=textprops,
                color=color,
            )
            ax.text(
                -0.75,
                i,
                source.report_fa_symbol,
                horizontalalignment="center",
                verticalalignment="center",
                fontproperties=fontawesome,
                color=color,
            )

        for assembly in assemblies_list:
            xa, ya = positions[assembly.id][0]
            if assembly.assembly_plan is not None:
                for _segment, child in assembly.assembly_plan.items():
                    xc, yc = positions[child.id][1]
                    ax.plot([xc, xa], [yc, ya], lw=1 * scale, color=color)

        ax.set_ylim((-1, len(assemblies_list) + 1))
        ax.set_xlim(-2, deadline)
        ax.set_frame_on(False)
        ax.yaxis.set_visible(False)
        ticks = list(np.arange(0, deadline, deadline / 6).astype(int))

        ax.xaxis.set_ticks(ticks)
        ax.xaxis.set_ticks_position("bottom")
        ax.xaxis.set_tick_params(which="major", length=5 * scale)
        ax.xaxis.set_ticklabels(ticks, fontproperties=textprops)
        ax.set_xlabel("TIME", fontproperties=textprops)
        ax.set_frame_on(False)
        return ax
