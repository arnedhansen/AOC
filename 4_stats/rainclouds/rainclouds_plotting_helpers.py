def _xpos_for_category_from_ticks(ax, category):
    """
    Infer the x position for a categorical label from current tick labels.
    This is a fallback when no explicit mapping is provided.
    """
    ticks = ax.get_xticks()
    labels = [t.get_text() for t in ax.get_xticklabels()]
    lab2x = {lab: x for lab, x in zip(labels, ticks)}
    if category not in lab2x:
        raise KeyError(f"Category '{category}' not found in axis tick labels: {labels}")
    return lab2x[category]


def add_stat_brackets(ax, xcats, comparisons, y_positions, labels,
                      bracket_height=0.02, lw=1.5, text_offset=0.01, fontsize=12,
                      xmap=None):
    """
    Draw significance brackets between category pairs on a categorical x-axis.

    Parameters
    ----------
    ax : matplotlib Axes
    xcats : sequence of category labels (for reference/order; not strictly required)
    comparisons : list of (cat1, cat2) tuples
    y_positions : list of float (data coords) for each comparison’s bracket baseline
    labels : list of strings (e.g., '*', 'n.s.')
    bracket_height : fraction of y-range used as bracket height
    lw : line width
    text_offset : fraction of y-range above bracket for the label
    fontsize : text size
    xmap : dict or None
        Optional explicit mapping {category_label: x_position}. If None, the function
        infers x from current xticks/xticklabels on the axes.
    """
    y0, y1 = ax.get_ylim()
    yr = (y1 - y0)
    h_px = bracket_height * yr
    tofs = text_offset * yr

    def _xpos_for_category(cat):
        if xmap is not None:
            if cat not in xmap:
                raise KeyError(f"Category '{cat}' not found in xmap keys: {list(xmap.keys())}")
            return xmap[cat]
        return _xpos_for_category_from_ticks(ax, cat)

    for (g1, g2), y, lab in zip(comparisons, y_positions, labels):
        x1 = _xpos_for_category(g1)
        x2 = _xpos_for_category(g2)
        if x1 > x2:
            x1, x2 = x2, x1

        ax.plot([x1, x1, x2, x2],
                [y,  y + h_px, y + h_px, y],
                linewidth=lw, color="black", clip_on=False)

        ax.text((x1 + x2) / 2.0, y + h_px + tofs, lab,
                ha="center", va="bottom", fontsize=fontsize)