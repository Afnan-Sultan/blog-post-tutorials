import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def _annotate_text(ax, text, loc=(0.5, 0.95), fontsize=10):
    ax.text(loc[0], loc[1], text,
            transform=ax.transAxes,
            fontsize=fontsize,
            ha='center',
            va='top')


def plot_distribution(series: pd.Series,
                              ax: plt.Axes,
                              x_label: str = '',
                              annotate_count: bool = True,
                              is_discrete: bool = False):
    ds_values = series.dropna()
    bin_edges = np.histogram_bin_edges(ds_values, bins='auto')
    sns.histplot(ds_values, bins=bin_edges, kde=not is_discrete, ax=ax)

    if not is_discrete:
        ax_box = ax.twinx()
        ax_box.set_ylim(0, 1)
        ax_box.set_xlim(ax.get_xlim())
        y_pos = 0.2
        ax_box.boxplot(ds_values,
                       vert=False,
                       positions=[y_pos],
                       widths=0.1,
                       patch_artist=True,
                       boxprops=dict(facecolor='none', alpha=0.7),
                       medianprops=dict(color='black', linewidth=1))
        ax_box.plot(np.mean(ds_values), y_pos, '^', color='green', ms=5)
        ax_box.set(yticks=[], ylabel='')

    if annotate_count:
        sample_count = len(ds_values)
        _annotate_text(ax,
                       text=f"# Molecules = {sample_count}",
                       fontsize=14,
                       loc=(0.2, 0.95))

    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.set_title(series.name, fontsize=18)
    # Increase tick label font sizes
    ax.tick_params(axis='both', labelsize=12)
