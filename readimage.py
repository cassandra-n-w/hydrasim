# This file is part of Hydrasim, a small piece of software intended to 
# simulate observations of TW Hydra.

# Hydrasim is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# Hydrasim is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with Hydrasim. If not, see <https://www.gnu.org/licenses/>. 


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
