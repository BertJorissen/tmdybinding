## TMDybinding
<img src="https://github.com/BertJorissen/tmdybinding/blob/master/docs/assets/images/logo.png?raw=true" width="100">

**TMDybinding** is an open-source Python module for [Pybinding].
This package contains some `pb.Lattice()`s to create tight-binding models for TMDs.

## Installation
Installation is simple, just run
```
pip install tmdybinding
```

## Usage
Below is a small example of how to use the code for the 5-band model from Liu/Wu:
```python
import pybinding as pb
import tmdybinding as td
from tmdybinding.sg_parameters import _params_sg_wu as params

# get the parameters
lat = td.TmdNN256Meo(params=params).lattice()
model = pb.Model(lat, pb.translational_symmetry())
bz = lat.brillouin_zone()
bands = pb.solver.lapack(model).calc_bands(
    bz[3] * 0, bz[3], (bz[3] + bz[4]) / 2, bz[3] * 0
)
bands.plot(point_labels=[r"$\Gamma$", "$K$", "$M$", r"$\Gamma$"])
```
<div>
  <figure>
    <img src="assets/images/bs_tmd.png" style="width: 60em; max-width: 50%; display: inline-block;"/>
    <figcaption>The band structure for the TMD.</figcaption>
  </figure>
</div>

Further documentation is coming soon.
Made by [Bert Jorissen].

[Bert Jorissen]: https://bertjorissen.be