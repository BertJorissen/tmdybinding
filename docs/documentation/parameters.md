## Possible parameters

In the subpackage `tmdybinding.parameters` you can find the parameters for the
Symmetry Group (SG) and Slater-Koster (SK) models, respectively. These parameters are dictionaries with the
index the name of the different TMDs considered and have as item the parameters in the form of a
[ParametersList](api.md#tmdybinding.ParametersList) object.

Here is an overview of the possible choices for the parameters and the best associated class to use for the parameter
from literature:

| Name                 | Type | Materials                                                   | Class               | Ref                      |
|----------------------|------|-------------------------------------------------------------|---------------------|--------------------------|
| `tmdy.liu2`          | SG   | `"MoS2"`, `"MoSe2"`, `"MoTe2"`, `"WS2"`, `"WSe2"`, `"WTe2"` | [TmdNN2Me]          | [G-B. Liu, 2013]         |
| `tmdy.liu6`          | SG   | `"MoS2"`, `"MoSe2"`, `"MoTe2"`, `"WS2"`, `"WSe2"`, `"WTe2"` | [TmdNN256Me]        | [G-B. Liu, 2013]         |
| `tmdy.wu`            | SG   | `"MoS2"`                                                    | [TmdNN256Meo]       | [F. Wu, 2015]            |
| `tmdy.fang`          | SG   | `"MoS2"`                                                    | [TmdNN123MeoXeo]    | [S. Fang, 2015]          |
| `tmdy.jorissen`      | SG   | `"MoS2"`                                                    | [TmdNN12MeXe]       | [B. Jorissen, 2024]      |
| `tmdy.all`           | SG   | `"MoS2"`                                                    | [TmdNN123456MeoXeo] | [B. Jorissen, 2024]      |
| `tmdy.rostami`       | SK   | `"MoS2"`                                                    | [TmdNN12MeXe]       | [H. Rostami, 2015]       |
| `tmdy.dias`          | SK   | `"MoS2"`, `"MoSe2"`, `"MoTe2"`, `"WS2"`, `"WSe2"`           | [TmdNN125MeoXeo]    | [A. Dias, 2018]          |
| `tmdy.cappelluti`    | SK   | `"MoS2"`                                                    | [TmdNN12MeoXeo]     | [E. Cappelluti, 2013]    |
| `tmdy.roldan`        | SK   | `"MoS2"`                                                    | [TmdNN12MeoXeo]     | [R. Roldan, 2014]        |
| `tmdy.ridolfi`       | SK   | `"MoS2": {"normal", "vb", "minimal"}`                       | [TmdNN12MeoXeo]     | [E. Ridolfi, 2015]       |
| `tmdy.venkateswarlu` | SK   | `"MoS2"`                                                    | [TmdNN12MeoXeo]     | [S. Venkateswarlu, 2020] |
| `tmdy.silva_guillen` | SK   | `"MoS2"`, `"MoSe2"`, `"WS2"`, `"WSe2"`                      | [TmdNN12MeoXeo]     | [J. Silva-Guillen, 2016] |
| `tmdy.pearce`        | SK   | `"MoS2"`                                                    | [TmdNN12MeoXeo]     | [A. Pearce, 2016]        |
| `tmdy.bieniek`       | SK   | `"MoS2"`                                                    | [TmdNN12MeoXeo]     | [M. Bieniek, 2018]       |
| `tmdy.abdi`          | SK   | `"MoS2"`,  `"MoSe2"`, `"WS2"`                               | [TmdNN12MeoXeo]     | [M. Abdi, 2020]          |
| `tmdy.peng6`         | SK   | `"MoS2"`, `"MoSe2"`, `"WS2"`, `"WSe2"`                      | [TmdNN12MeXe]       | [Z. Peng, 2024]          |
| `tmdy.peng11`        | SK   | `"MoS2"`, `"MoSe2"`, `"WS2"`, `"WSe2"`                      | [TmdNN125MeoXeo]    | [Z. Peng, 2024]          |

!!! warning "**Selecting other classes**"

    The fourth column from the table above is the suggested model to be used with the parameters. However, you can use
    any of the classes to build the matrices associated to the hoppings of that model. The only difference is that some
    parameters that are not orginally defined for that model will be set to zero, giving unexpected results.


!!! tip "**Starting from Slater-Koster parameters**"

    The TMDybinding package automatically converts the Slater-Koster parameters to the Symmetry Group parameters.
    If you save your variables in the [`SKParametersList`](api.md#tmdybinding.SKParametersList) object, you can obtain the
    Symmetry Group parameters by calling the appropriate parameter to this class.

    You also have the option to specify the angle between the plane of the Metal atoms and the first nearest neighbour
    hopping direction to the Chalcogenide atoms.


[comment]: # (The following are the links to the classes)
[TmdNN2Me]: api.md#tmdybinding.TmdNN2Me
[TmdNN256Me]: api.md#tmdybinding.TmdNN256Me
[TmdNN256Meo]: api.md#tmdybinding.TmdNN256Meo
[TmdNN123MeoXeo]: api.md/#tmdybinding.TmdNN123MeoXeo
[TmdNN12MeXe]: api.md#tmdybinding.TmdNN12MeXe
[TmdNN123456MeoXeo]: api.md#tmdybinding.TmdNN123456MeoXeo
[TmdNN125MeoXeo]: api.md#tmdybinding.TmdNN125MeoXeo
[TmdNN12MeoXeo]: api.md#tmdybinding.TmdNN12MeoXeo

[comment]: # (The following are the links to the publications)
[G-B. Liu, 2013]: https://doi.org/10.1103/PhysRevB.88.085433
[F. Wu, 2015]: https://doi.org/10.1103/PhysRevB.91.075310
[S. Fang, 2015]: https://doi.org/10.1103/PhysRevB.92.205108
[B. Jorissen, 2024]: https://doi.org/10.21468/SciPostPhysCore.7.1.004
[H. Rostami, 2015]: https://doi.org/10.1103/PhysRevB.92.195402
[A. Dias, 2018]: https://doi.org/10.1103/PhysRevB.98.075202
[E. Cappelluti, 2013]: https://doi.org/10.1103/PhysRevB.88.075409
[R. Roldan, 2014]: https://doi.org/10.1088/2053-1583/1/3/034003
[E. Ridolfi, 2015]: https://doi.org/10.1088/0953-8984/27/36/365501
[S. Venkateswarlu, 2020]: https://doi.org/10.1103/PhysRevB.102.081103
[J. Silva-Guillen, 2016]: https://doi.org/10.3390/app6100284
[A. Pearce, 2016]: https://doi.org/10.1103/PhysRevB.94.155416
[M. Bieniek, 2018]: https://doi.org/10.1103/PhysRevB.97.085153
[M. Abdi, 2020]: https://doi.org/10.1149/2162-8777/abb28b
[Z. Peng, 2024]: https://doi.org/10.1103/PhysRevB.109.245412
