import os
import numpy as np
import matplotlib.pyplot as plt

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def _xpos_for_category(ax, cat):
    # map categorical tick labels to x positions
    labels = [tick.get_text() for tick in ax.get_xticklabels()]
    if cat not in labels:
        # force draw to populate tick labels if needed
        plt.draw()
        labels = [tick.get_text() for tick in ax.get_xticklabels()]
    try:
        idx = labels.index(cat)
    except ValueError:
        # fall back to numeric position by scanning artists (robust enough)
        idx = 0
    return idx

def add_stat_brackets(ax, xcats, comparisons, y_positions, labels,
                      bracket_height=0.02, lw=1.5, text_offset=0.01, fontsize=12):
    # Draws significance brackets between category pairs on a categorical x-axis
    for (g1, g2), y, lab in zip(comparisons, y_positions, labels):
        x1 = _xpos_for_category(ax, g1)
        x2 = _xpos_for_category(ax, g2)
        if x1 > x2:
            x1, x2 = x2, x1

        # width of bracket
        h = bracket_height * (ax.get_ylim()[1] - ax.get_ylim()[0])

        # bracket lines
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], linewidth=lw, color="black", clip_on=False)

        # text label (asterisks)
        ax.text((x1 + x2) / 2, y + h + text_offset * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                lab, ha="center", va="bottom", fontsize=fontsize)
