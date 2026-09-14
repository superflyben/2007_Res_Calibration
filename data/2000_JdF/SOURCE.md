# 2000 Juan de Fuca calibration data

Five isohaline runs from the 2000 Juan de Fuca bench calibration, 20 points total.
Columns are unlabelled in the original files:

1. NaCl concentration, **mol/kg** (note: mol, not mmol -- multiply by 1000 to match
   the 2007 data)
2. Temperature, deg C
3. Conductivity, 1/V

Pressure is 220 bar, except the 382 deg C block which was run at 300 bar.

Source: the calibration behind

> Larson, B. I., Olson, E. J. and Lilley, M. D. (2007), "In situ measurement of
> dissolved chloride in high temperature hydrothermal fluids", *Geochimica et
> Cosmochimica Acta*. doi:10.1016/j.gca.2007.02.013

Copied verbatim from `Udub/Modeling/3rdTSfit/`, where they are consumed by
`TSfit3.m` via `readlarson.m`. `TSfit3.m` builds its surface by spline-interpolating
each isohaline along temperature (`interpmmlarson.m`, 10000 points) and then gridding
with `gridfit` at 200x200. This project deliberately does **not** reuse that surface:
both calibrations are rebuilt from raw points with the same algorithm so that the
comparison between them reflects the data rather than the fitting method.
