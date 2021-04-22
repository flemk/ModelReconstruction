import numpy as np
import stanpy as sp

print('loading series')
time_series = np.load('VdP_series_xv.npy')
print('loading finished')

analysis = sp.StochasticAnalysis(time_series)
analysis.analyze(1)

np.save('VdP_coefficients_analysis.npy', analysis)
