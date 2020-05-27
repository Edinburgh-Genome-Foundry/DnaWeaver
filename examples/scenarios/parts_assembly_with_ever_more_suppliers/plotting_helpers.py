import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from dnaweaver.AssemblyPlanReport.mixins.PlotsMixin.plot_supply_network import (
    plot_supply_network,
)


def matplotlib_gridspecs_from_array(arr):
    arr = np.array(arr)
    h, w = arr.shape
    values_in_arr = sorted(set(arr.flatten()))
    result = []
    grid = GridSpec(h, w)
    for v in values_in_arr:
        xx, yy = ((arr - v) == 0).nonzero()
        result.append(grid[xx.min() : xx.max() + 1, yy.min() : yy.max() + 1])
    return result


def matplotlib_axes_from_gridspec_array(arr, figsize=None):
    """Returned axes layed out as indicated in the array

    Example:
    --------
    >>> # Returns 3 axes layed out as indicated by the array
    >>> fig, axes = matplotlib_axes_from_gridspec_array([
    >>>     [1, 1, 3],
    >>>     [2, 2, 3],
    >>>     [2, 2, 3],
    >>> ])
    """
    fig = plt.figure(figsize=figsize)
    gridspecs = matplotlib_gridspecs_from_array(arr)
    axes = []
    for gridspec in gridspecs:
        axes.append(fig.add_subplot(gridspec))
    return fig, axes


def plot_quote(quote, figsize=(4, 5), ylim=None, axes=None):
    """Plot a quote (supply network and assembly plan)"""

    report = quote.to_assembly_plan_report()
    if axes is None:
        fig, axes = matplotlib_axes_from_gridspec_array(
            [[1], [1], [2]], figsize=figsize
        )
    else:
        fig = axes[0].figure
    plot_supply_network(report, ax=axes[0], margin=0.2)
    report.plot_assembly_blocks(
        parts_offset=0.1, plot_top_assembly=False, ax=axes[1], legend=True
    )
    if ylim is not None:
        axes[1].set_ylim(*ylim)
    return fig, axes
