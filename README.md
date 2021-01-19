# ModelReconstruction
Analysis of time series from stochastic processes.

## ```stanpy.py``` - Stochystical Analysis Python
> Based on "Analysis of time series from stochastic processes" (by J.Gradisek, S.Siegert, R.Friedrich, I.Grabec)

> ```StochasticAnalysis.ipynb``` shows derivation of the code and some examples.

This module provides a class to determine drift- and diffusion-coefficients of n-dimensional time series by using their statistical definition.

### Usage
```python
import stanpy as sp

time_series = [[1, 2, ...], [1, 2, ...]] # your time seres you want to analyze

analysis = sp.StochasticAnalysis(time_series)
analysis.analyze()

# drift and diffusion coefficients are now stored in:
analysis.drift()
analysis.diffusion()

# in the 2d case you can visualize them builtin:
analysis.visualize_2d()

# and you can reconstruct your series with choosen initial values:
r = analysis.reconstruct()

# by converting your coefficients into a FPE you might gain more insight:
f = analysis.solve_fpe()
```

## ```FPE.ipynb``` - solution of the Fokker-Planck-Equation
The FPE of two cases (Harmonic- and Van-der-Pol-oscilaltor) is being solved by euler integration.

## ```casesBW.ipynb``` - analyzing covid-19 data of Baden-Württemberg, Germany
Title says it all.