The TMDybinidng package consists of several base lattices, objects to store the parameters for these lattices,
functions to build the matrices associated to the hoppings of that model and classes that help user in
defining new lattices with (rotational) symmetric aspects for orbitals of non-zero ml-quantum number.

## Possible parameters

In the files `tmdybinding/sg_parameters.py` and `tmdybinding/sk_parameters.py` you can find the parameters for the
Symmetry Group (SG) and Slater-Koster (SK) models, respectively. These parameters are dictionaries with the
index the name of the different TMDs considered and have as item the parameters in the form of a
[ParametersList](#tmdybinding.ParametersList) object.

Here is an overview of the possible choices for the parameters and the best associated class to use for the parameter
from literature:

| Name                 | Type | Materials                                                   | Class                                               | Ref                                                                                  |
|----------------------|------|-------------------------------------------------------------|-----------------------------------------------------|--------------------------------------------------------------------------------------|
| `tmdy.liu2`          | SG   | `"MoS2"`, `"MoSe2"`, `"MoTe2"`, `"WS2"`, `"WSe2"`, `"WTe2"` | [TmdNN2Me](#tmdybinding.TmdNN2Me)                   | [G-B. Liu, 2013](https://link.aps.org/doi/10.1103/PhysRevB.88.085433)                |
| `tmdy.liu6`          | SG   | `"MoS2"`, `"MoSe2"`, `"MoTe2"`, `"WS2"`, `"WSe2"`, `"WTe2"` | [TmdNN256Me](#tmdybinding.TmdNN256Me)               | [G-B. Liu, 2013](https://link.aps.org/doi/10.1103/PhysRevB.88.085433)                |
| `tmdy.wu`            | SG   | `"MoS2"`                                                    | [TmdNN256Meo](#tmdybinding.TmdNN256Meo)             | [F. Wu, 2015](https://link.aps.org/doi/10.1103/PhysRevB.91.075310)                   |
| `tmdy.fang`          | SG   | `"MoS2"`                                                    | [TmdNN123MeoXeo](#tmdybinding.TmdNN123MeoXeo)       | [S. Fang, 2015](https://link.aps.org/doi/10.1103/PhysRevB.92.205108)                 |
| `tmdy.jorissen`      | SG   | `"MoS2"`                                                    | [TmdNN12MeXe](#tmdybinding.TmdNN12MeXe)             | [B. Jorissen, 2024](https://scipost.org/SciPostPhysCore.7.1.004)                     |
| `tmdy.all`           | SG   | `"MoS2"`                                                    | [TmdNN123456MeoXeo](#tmdybinding.TmdNN123456MeoXeo) | [B. Jorissen, 2024](https://scipost.org/SciPostPhysCore.7.1.004)                     |
| `tmdy.rostami`       | SK   | `"MoS2"`                                                    | [TmdNN12MeXe](#tmdybinding.TmdNN12MeXe)             | [H. Rostami, 2015](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.195402) |
| `tmdy.dias`          | SK   | `"MoS2"`, `"MoSe2"`, `"MoTe2"`, `"WS2"`, `"WSe2"`           | [TmdNN125MeoXeo](#tmdybinding.TmdNN125MeoXeo)       |                                                                                      |
| `tmdy.cappelluti`    | SK   | `"MoS2"`                                                    | [TmdNN12MeoXeo](#tmdybinding.TmdNN12MeoXeo)         |                                                                                      |
| `tmdy.roldan`        | SK   | `"MoS2"`                                                    | [TmdNN12MeoXeo](#tmdybinding.TmdNN12MeoXeo)         |                                                                                      |
| `tmdy.ridolfi`       | SK   | `"MoS2": {"normal", "vb", "minimal"}`                       | [TmdNN12MeoXeo](#tmdybinding.TmdNN12MeoXeo)         |                                                                                      |
| `tmdy.venkateswarlu` | SK   | `"MoS2"`                                                    | [TmdNN12MeoXeo](#tmdybinding.TmdNN12MeoXeo)         |                                                                                      |
| `tmdy.silva_guillen` | SK   | `"MoS2"`, `"MoSe2"`, `"WS2"`, `"WSe2"`                      | [TmdNN12MeoXeo](#tmdybinding.TmdNN12MeoXeo)         |                                                                                      |
| `tmdy.bieniek`       | SK   | `"MoS2"`                                                    | [TmdNN12MeoXeo](#tmdybinding.TmdNN12MeoXeo)         |                                                                                      |
| `tmdy.abdi`          | SK   | `"MoS2"`,  `"MoSe2"`, `"WS2"`                               | [TmdNN12MeoXeo](#tmdybinding.TmdNN12MeoXeo)         |                                                                                      |

::: tmdybinding
