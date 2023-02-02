#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.cm as cm
import numpy as np
from matplotlib.colors import ListedColormap


class MatplotlibConfig:
    """
    A class with color and structure specifics
    """

    top = cm.get_cmap("Blues_r", 128)
    bottom = cm.get_cmap("Oranges", 128)
    newcolors = np.vstack((top(np.linspace(0, 1, 128)), bottom(np.linspace(0, 1, 128))))
    newcmp = ListedColormap(newcolors, name="OrangeBlue")


color = {
    "Core proteome genes": MatplotlibConfig.newcmp(0.85),
    "Other proteins": MatplotlibConfig.newcmp(0.15),
}
