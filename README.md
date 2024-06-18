import tmdybinding as td

# make the lattice object
lat = td.TmdNN256Me(params=td.liu6["MoS2"]).lattice()

# make a figure with the lattice
lat.plot()# TMDybinding
<img src="https://github.com/BertJorissen/tmdybinding/blob/master/docs/assets/images/logo.png?raw=true" width="100">

**TMDybinding** is an open-source Python module for [Pybinding].

## Capabilities

* generate accurate lattices for most used TB models for TMDs
* easily add hopping
* automatic generation of lattice cell with 90-degree corners
* automatic inclusion of SOC terms
* automatic generation of hopping terms taking into account the symmetry of the orbital
* label all the different orbital contributions
* fully compatible with [pybinding] and [KITE]

## Installation

That's easy, just run
``` 
pip install tmdybinding
```


[pybinding]: https://pybinding.site
[KITE]: https://quantum-kite.com