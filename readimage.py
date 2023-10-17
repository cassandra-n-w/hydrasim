# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
from matplotlib import cm
from radmc3dPy import *
im=image.readImage()
image.plotImage(im,vmax=3e-3,au=True,cmap=cm.gist_heat)
