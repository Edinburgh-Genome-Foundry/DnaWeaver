import tqdm
import time
import matplotlib.pyplot as plt
from matplotlib import rc, font_manager, rcParams
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy
from generate_supply_network import generate_supply_network


# GENERATE THE DATA


def generate_quote(a_star_factor):
    t0 = time.time()
    station = generate_supply_network(a_star_factor)
    station.prepare_network_on_sequence(sequence)
    quote = station.get_quote(sequence, with_assembly_plan=True)
    t1 = time.time()
    return quote, t1 - t0


with open("sequence.txt", "r") as f:
    sequence = f.read()

a_star_factors = numpy.logspace(-3, 0.5, 20)
quotes = [
    (a_start_factor, [generate_quote(a_start_factor) for i in range(5)])
    for a_start_factor in tqdm.tqdm(a_star_factors)
]
average_times = [numpy.mean([t for _, t in data]) for a, data in quotes]
prices = [data[0][0].price for a, data in quotes]


# PLOT THE PRICE AND TIME CURVES

fig, ax = plt.subplots(1, figsize=(7, 3))
ax2 = ax.twinx()

dark_blue = "#297695"
dark_red = "#e75c5c"
ax.plot(a_star_factors, average_times, marker="o", c=dark_blue, lw=2)
ax2.plot(a_star_factors, prices, marker="o", c=dark_red, lw=2)
ax.set_xscale("log")
ax2.set_ylim(bottom=0, top=10000)
ax.axvline(x=0.15, ls=":", c="k")
ax.axvline(x=0.15, ls=":", c="k")
ax.set_yticks([0.8, 16])
ax.set_yticklabels(["0.8s", "16s"], fontdict=dict(size=14, color=dark_blue))
ax.set_xticks([0.01, 0.1, 1])
ax.set_xticklabels(["1c/bp", "10c/bp", "£1/bp"], fontdict=dict(size=14))

yticks = [quotes[0][1][0][0].price, quotes[-1][1][0][0].price]
ax2.set_yticks(yticks)
ax2.set_yticklabels(
    ["£" + "%d" % c for c in yticks], fontdict=dict(size=14, color=dark_red)
)
ax.spines["top"].set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.set_xlabel("A* factor")
fig.tight_layout()
fig.savefig("a_star_factor_effect.svg")
fig.savefig("a_star_factor_effect.pdf")

# PLOT THE BLOCKS

fig, axes = plt.subplots(1, 3, figsize=(12, 3))
selected_quotes = [quotes[i] for i in (0, 12, 19)]
for ax, (factor, quotes_list) in zip(axes, selected_quotes):
    quote = quotes_list[0][0]
    mean_time = numpy.mean([t for (q, t) in quotes_list])
    report = quote.to_assembly_plan_report()
    report.plot_assembly_blocks(
        ax=ax,
        parts_offset=0.1,
        plot_top_assembly=False,
        legend=(ax == axes[0]),
        legend_offset=-0.05,
    )
    ax.set_title(
        "A* factor %.02f - £%d\nComputed in %.1fs"
        % (factor, quote.price, mean_time)
    )
    if ax.legend_ is not None and (ax != axes[0]):
        ax.legend_.remove()

fig.savefig("a_star_factor_blocks.svg")
