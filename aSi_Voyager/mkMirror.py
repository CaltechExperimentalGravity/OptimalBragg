#!/usr/bin/env python

import sys
sys.path.append('../../pygwinc/')


sys.path.append('../generic/pythonAddOns/')
from optimUtils import *

ifo = gwinc.load_ifo('aSiModel.m')

voy = gwinc.load_ifo('Voyager')

Ls = np.array(voy.Optics.ETM.CoatLayerOpticalThickness)

gam = brownianProxy(ifo)

getMirrorCost('params.yml', Ls, ifo, gam, verbose=True)
