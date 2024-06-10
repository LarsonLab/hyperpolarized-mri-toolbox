"""
This module was written to store potentially useful visualization tools for
hp data (aka high dimensional image data)
author: anna bennett
date: 05.30.2024
"""

import matplotlib.pyplot as plt
import numpy as np


def multiframeimshow(images, scale, dims, customCMAP=None, mask=None):
    fig, axs = plt.subplots(
        dims[0], dims[1], facecolor="w", edgecolor="k", tight_layout=True
    )

    if customCMAP is None:
        colors = "viridis"
    else:
        colors = customCMAP

    if mask is not None:
        images[mask == 0] = np.nan

    for ax, i in zip(axs.ravel(), range(images.shape[2])):
        ax.set_xlim(0, images.shape[0] - 1)
        ax.set_ylim(0, images.shape[1] - 1)
        ax.imshow(
            images[:, :, i],
            vmin=scale[0],
            vmax=scale[1],
            origin="upper",
            cmap=colors,
            extent=[0, images.shape[0], 0, images.shape[1]],
            aspect=1,
        )

        y_ax = ax.axes.get_yaxis()
        y_ax.set_visible(False)
        x_ax = ax.axes.get_xaxis()
        x_ax.set_visible(False)

    return fig
