'''
Stochastic Analysis of time series

Franz Ludwig Kostelezky, 2020
'''

import numpy as np
import matplotlib.pyplot as plt

class StochasticAnalysis:
    def __init__(self, series, drift=None, diffusion=None, reconstruction=None):
        self._series = series
        self._drift = drift
        self._diffusion = diffusion
        self._recontruction = reconstruction

        self._fpe = None
        self._dimension = len(series)

    def D_1(self, series, dt=1, bins=250, tau=1, transform=None):
        ''' Retrieving n-dimensional Drift-Coefficient

        Parameters:
            - (array) series: array of n arrays which represent time series. n-dimensinal.
            - (float) dt: (time) difference between the values of series.
            - (int) bins: Number of bins. Defines the accuracy.
            - (int) tau: Number of timesteps to derivate further.
            - (lambda) transform: transform function to project series value to mesh.
                default: transform = lambda x, d, b: int((x + (d/2)) * np.floor(b / d)) - 1

        Returns:
            (array) Array of n-dim arrays, where the arrays represent the mean change of the i-th variable.
        '''
        if self._drift is not None:
            print('Drift is already defined. Will be overwritten by newer value.')

        dimension = len(series)

        # checking if all series have same size
        for i in range(dimension):
            if len(series[i]) != len(series[0]):
                raise Exception('Not all series have the same length')

        d = [np.max(el) - np.min(el) for el in series] # offsets
        l = [np.zeros(bins) for _ in range(dimension)] # n-dimension array

        a_grid = np.meshgrid(*l) # mesh to store changes
        b_grid = np.meshgrid(*l) # mesh to count occurences

        if transform is None:
            transform = lambda x, d, b: int((x + (d/2)) * np.floor(b / d)) - 1

        for i in range(len(series[0][:-tau])):
            # 1. transform series value to index value of grids
            c = [transform(series[j][i], d[j], bins) for j in range(dimension)]
            c = tuple(c)
            # 2. summate changes of the series and write to mesh
            for j in range(dimension):
                a_grid[j][c] += series[j][i + tau] - series[j][i]
            # 3. increment number of visits
            for j in range(dimension):
                b_grid[j][c] += 1

        # now calculate mean changes
        def calculate_mean_change_recursive(s, s_, argument=1):
            if type(s) == np.ndarray:
                for i in range(len(s)):
                    s[i] = calculate_mean_change_recursive(s[i], s_[i], argument=argument)
            else:
                result = 0 if s_ == 0 else s / s_
                return argument * result
            return s
        for j in range(dimension):
            a_grid[j] = calculate_mean_change_recursive(a_grid[j], b_grid[j], argument=(1 / (tau * dt)))

        # save in class for later usage
        self._drift = a_grid

        return a_grid

    def drift(self):
        ''' Returns calculated Drift-coefficitent.
        '''
        return self._drift

    def D_2(self, series, dt=1, bins=250, tau=1, transform=None):
        ''' Retrieving n-dimensional Drift-Coefficient

        Parameters:
            - (array) series: array of n arrays which represent time series. n-dimensinal.
            - (float) dt: (time) difference between the values of series.
            - (int) bins: Number of bins. Defines the accuracy.
            - (int) tau: Number of timesteps to derivate further.
            - (lambda) transform: transform function to project series value to mesh.
                default: transform = lambda x, d, b: int((x + (d/2)) * np.floor(b / d)) - 1

        Returns:
            (array) Array of n-dim arrays, where the arrays represent the mean change of the i-th variable.
        '''
        if self._diffusion is not None:
            print('Diffusion is already defined. Will be overwritten by newer value.')

        dimension = len(series)

        # checking if all series have same size
        for i in range(dimension):
            if len(series[i]) != len(series[0]):
                raise Exception('Not all series have the same length')

        d = [np.max(el) - np.min(el) for el in series] # offsets
        l = [np.zeros(bins) for _ in range(dimension)] # n-dimension array

        a_grid = np.meshgrid(*l) # mesh to store changes
        b_grid = np.meshgrid(*l) # mesh to count occurences

        a_grid = [a_grid[0] for _ in range(dimension * dimension)]

        if transform is None:
            transform = lambda x, d, b: int((x + (d/2)) * np.floor(b / d)) - 1

        for i in range(len(series[0][:-tau])):
            # 1. transform series value to index value of grids
            c = [transform(series[j][i], d[j], bins) for j in range(dimension)]
            c = tuple(c)
            # 2. summate and multiply changes of the series and write to mesh
            for k in range(dimension):
                for j in range(dimension):
                    d_c = k * dimension + j
                    a_grid[d_c][c] += (series[j][i + tau] - series[j][i]) * (series[k][i + tau] - series[k][i])
            # 3. increment number of visits
            for j in range(dimension):
                b_grid[j][c] += 1

        # now calculate mean changes
        def calculate_mean_change_recursive(s, s_, argument=1):
            if type(s) == np.ndarray:
                for i in range(len(s)):
                    s[i] = calculate_mean_change_recursive(s[i], s_[i], argument=argument)
            else:
                result = 0 if s_ == 0 else s / s_
                return argument * result
            return s
        for j in range(dimension):
            a_grid[j] = calculate_mean_change_recursive(a_grid[j], b_grid[j], argument=(1 / (tau * dt)))

        # save in class for later usage
        self._diffusion = a_grid

        return a_grid

    def diffusion(self):
        ''' Returns calculated Diffusion-coefficitent.
        '''
        return self._diffusion

    def recontruct(self, ivp=None):
        ''' Reconstruct a time series with given initial values based on the precalculated coefficients

        Parameters:
            - (np.ndarray) ivp: Initial values in shape np.shape(series[0]). If None np.zeros(np.shape(series[0])) is choosen.

        Returns:
            - (np.ndarray) reconstructed: Calculated reconstructed timeseries.
        '''
        pass

    def solve_fpe(self):
        ''' Series is converted to a Fokker-Planck-equation and then solved.
        '''
        pass

    def visualize_2d(self, drift, diffusion):    
        ''' Plots precalculated drift- and diffusion-coefficients.

        Returns:
            (matplotlib.pyplot.figure) 2D tuple of figures.
        '''

        # abort if dimension is not 2 - not plotable
        if self._dimension != 2:
            return 0

        figs = []

        # add drift coefficient D_1 to plot
        fig = plt.figure()

        dimension = len(drift)
        for i in range(dimension):
            plt.subplot(1, dimension, i + 1)
            plt.imshow(drift[i], cmap='hot')
            plt.title('$D^{(1)}_{x_' + str(i) + '}$')
            plt.xlabel('$x_0$')
            plt.ylabel('$x_1$')

        figs.append(fig)

        # add diffusion coefficient D_2 to plot
        fig = plt.figure()

        a = lambda x, dim: int((((x - (x % dim)) / dim)) * dim)
        b = lambda x, dim: int(x - a(x, dim))

        dimension = len(diffusion)
        sqrtdim = int(np.sqrt(dimension))
        for i in range(dimension):
            plt.subplot(1, dimension, i + 1)
            plt.imshow(diffusion[i], cmap='hot')
            plt.title('$D^{(2)}_{x_{' + str(a(i, sqrtdim)) + '}x_{' + str(b(i, sqrtdim)) + '}}$')
            plt.xlabel('$x_0$')
            plt.ylabel('$x_1$')

        figs.append(fig)

        return figs

    def analyze(self, dt):
        ''' Calculate drift and diffusion coefficients.
        ... using only one command.
        '''
        _ = self.D_1(self._series, dt)
        _ = self.D_2(self._series, dt)

        return 0

class FPE:
    def __init__(self, drift, diffusion, ivp, dx):
        self._drift = drift
        self._diffusion = diffusion
        self._ivp = ivp
        self._dx = dx

    def upwind_scheme(self, x, n_order=1):
        ''' Returns the dimensional derivate in n_order - order upwind-scheme.
        Particularly programmed for probability field W in the Fokker-Planck-Equation.

        Parameters:
            - (np.ndarray) x: vector or field to derivate
            - (int) n_order: order of upwind-scheme. n-order in [1, 2, 3]

        Returns:
            - (np.ndarray) x_: vector derivate
        '''
        pass

    def solve(self, t_interval, dt):
        ''' Solve the Fokker-Planck-Equaion by euler-integration.
        '''
        pass

# 1. standard transform function
transform = lambda x, d, b: int((x + (d/2)) * np.floor(b / d)) - 1

# 2. standard axis
normaxis = lambda bins, dimension: np.meshgrid(*[np.arange(0, bins) for _ in range(dimension)])
