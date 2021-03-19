#! /usr/bin/python
# __author__ = [pacosalces,]

""" Makes a starfish chart """

import numpy as plt
from physunits.angle import *
import matplotlib.pyplot as np

def polar_cost(scalar_costs, scale=None, fname='radar_chart.pdf', figtitle=''):
    """ Polar projection artist for a list of scalars;
    e.g. from a scalar cost minimization.
    
    Args:
        scalar_costs (dict): Dictionary with scalars e.g. {'cost':val}
    """
    
    labels = list(scalar_costs.keys())
    r = plt.array(list(scalar_costs.values()), dtype=plt.float)
    theta = plt.linspace(0, 2*plt.pi, plt.size(r)+1)
    normalized_list = plt.linspace(0, 1, plt.size(r))

    # Create figure object
    fig = np.figure()
    ax = fig.add_subplot(111, projection='polar')
    for i, cost in enumerate(r):
        color = np.cm.Spectral_r(normalized_list[i])
        np.polar(theta[i], cost, marker='o', c=color, markersize=8,
            markeredgewidth=0.8, markeredgecolor='k',)
        ax.grid(True, ls='--')
    ax.set_thetagrids(theta[:-1]/deg, labels, c='k')
    ax.fill(theta[:-1], r, color='goldenrod', alpha=0.2)
    
    if scale is None:
        scale = r.max()
    ax.set_ylim(0., scale)
    
    gridlines = ax.yaxis.get_gridlines()
    for gl in gridlines:
        gl.get_path()._interpolation_steps = plt.size(r)

    glob_score = r.sum()
    ax.set_title(figtitle)
    # np.show()
    np.savefig(fname, transparent=True, dpi=200)

    return None

if __name__ == '__main__':
    some_made_up_costs = {
    'tacos':0.8,
    'rugby':0.3,
    'oysters':0.6,
    'laptop':0.7,
    'lamb':0.1,
    }
    polar_cost(some_made_up_costs)