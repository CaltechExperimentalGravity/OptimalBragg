# Python to calculate the cost function for a given coating

The functions are in `generic_local/optimUtils.py`. 
Example usage:
```python
>>> import sys
>>> from generic_local.optimUtils import *
>>> ifo = gwinc.load_ifo('Ta2O5Model.m')
>>> voy = gwinc.load_ifo('Voyager')
>>> Ls = np.array(voy.Optics.ETM.CoatLayerOpticalThickness)
>>> gam = brownianProxy(ifo)
>>> getMirrorCost('params.yml', Ls, ifo, gam, verbose=True)
```
The file `params.yml` specifies some properties for the cost func calculation (e.g. what terms to use, weights, number of layer pairs, etc..).

To do:
- [] Make a cleaner way of specifying which costs are computed (right now it's not completely modular)
- [] Merge dichroic enabling changes from optimUtils.py back into the `../generic/` directory.


