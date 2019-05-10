# Python to calculate the cost function for a given coating

The functions are in `Coatings/generic/pythonAddOns/optimUtils.py`. 
Example usage (I tested this against the MATLAB cost function and got the same cost):
```python
>>> from optimUtils import *
>>> ifo = gwinc.load_ifo('aSiModel.m')
>>> voy = gwinc.load_ifo('Voyager')
>>> Ls = np.array(voy.Optics.ETM.CoatLayerOpticalThickness)
>>> gam = brownianProxy(ifo)
>>> getMirrorCost('params.yml', Ls, ifo, gam, verbose=True)
```
The file `params.yml` specifies some properties for the cost func calculation (e.g. what terms to use, weights)

