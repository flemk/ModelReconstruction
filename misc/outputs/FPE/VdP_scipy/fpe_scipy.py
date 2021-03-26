import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from scipy.integrate import solve_ivp

from custom_utils import *

# initialize ivp and coefficients
bins = 50
interval = (-4, 4)
dx = (max(interval) - min(interval)) / bins

x_1, x_2 = np.meshgrid(np.linspace(min(interval), max(interval), bins), np.linspace(min(interval), max(interval), bins))
wxv = np.zeros((bins, bins))

# This choice of parameters follows Gadisek, et.al. PRE62, p. 3155 (2000)
mu = 1
epsilon = 2
omega = 1
dt = 0.00001
g = 1

d_x = x_2
d_v = mu * (epsilon - x_1 ** 2) * x_2 - (omega ** 2) * x_1

#wxv[75:125, 145:195] = 100 # "rechts"
#wxv[75:125, 75:125] = 100 # "mitte"
#wxv[75:125, 5:55] = 100 # "links"
wxv = np.exp(- 2 * x_1 ** 2 - 2 * x_2 ** 2)
#wxv[3:13, 12:25] = 100

wxv_ivp = wxv # this is our ivp

# the coefficients upwind scheme
d_x_pos = np.array(d_x)
d_x_pos[d_x_pos < 0] = 0 # all negative values of d_x are replaced by 0

d_x_neg = np.array(d_x)
d_x_neg[d_x_neg >= 0] = 0 # all positive values of d_x are replaced by 0 including 0

d_v_pos = np.array(d_v)
d_v_pos[d_v_pos < 0] = 0 # all negative values of d_x are replaced by 0

d_v_neg = np.array(d_v)
d_v_neg[d_v_neg >= 0] = 0 # all positive values of d_x are replaced by 0 including 0

# now for the actual solving.
# This is our function
# y equals the 1d version of wxv:
y = wxv.flatten()

def fokker_planck_system(t, y):
    if 1 or (t % .1 == 0):
        print('currently at t:%f' % (t))

    # calculate wxv "field" from 1d flattend array y
    wxv = y
    wxv = wxv.reshape((bins, bins))

    # extended bounds
    # used so that derivated use only valid values
    wxv_extended = np.empty((bins + 4, bins + 4))
    wxv_extended[2:-2, 2:-2] = wxv
    # apply neumann boundary conditions
    # with extended bounds
    wxv_extended[1, :] = (4. * wxv_extended[2, :] - wxv_extended[3, :]) / 3.
    wxv_extended[-2, :] = (4. * wxv_extended[-3, :] - wxv_extended[-4, :]) / 3.
    wxv_extended[:, -2] = (4. * wxv_extended[:, -3] - wxv_extended[:, -4]) / 3.
    wxv_extended[:, 1] = (4. * wxv_extended[:, 2] - wxv_extended[:, 3]) / 3.
    # extended bounds
    # propagate the boundary conditions to bounds
    wxv_extended[0, :] = wxv_extended[1, :]
    wxv_extended[-1, :] = wxv_extended[-2, :]
    wxv_extended[:, -1] = wxv_extended[:, -2]
    wxv_extended[:, 0] = wxv_extended[:, 1]
    # at the corners
    # also propagated to outer bound
    wxv_extended[1, 1] = 4. * wxv_extended[2, 2] - wxv_extended[3, 3]
    wxv_extended[0, 0] = wxv_extended[1, 1]
    wxv_extended[-2, 1] = 4. * wxv_extended[-3, 2] - wxv_extended[-4, 3]
    wxv_extended[-1, 0] = wxv_extended[-1, 1]
    wxv_extended[1, -2] = 4. * wxv_extended[2, -3] - wxv_extended[3, -4]
    wxv_extended[0, -1] = wxv_extended[1, -2]
    wxv_extended[-2, -2] = 4. * wxv_extended[-3, -3] - wxv_extended[-4, -4]
    wxv_extended[-1, -1] = wxv_extended[-2, -2]

    # receive first order upwind derivates
    el = first_order_upwind(wxv_extended, dx)
    # remove extended bounds
    el = np.array(el)
    el_ = []
    for i in range(len(el)):
        el_.append(el[i][2:-2, 2:-2])
    dwdx_pos, dwdx_neg, dwdv_pos, dwdv_neg = el_

    # the second derivate of diffusion coefficient is constantly
    d2wdv2 = (np.roll(wxv_extended, shift=(0, -1), axis=(1, 0)) \
              - 2 * wxv_extended + np.roll(wxv_extended, shift=(0, 1), axis=(1, 0))) / (dx ** 2)
    # remove extended bounds
    d2wdv2 = d2wdv2[2:-2, 2:-2]

    # now apply for each point in grid y
    y_ = np.zeros(np.shape(wxv))
    y_ = (+ (epsilon - x_1 ** 2) * wxv + (g / 2) * d2wdv2 \
          - d_x_pos * dwdx_pos - d_x_neg * dwdx_neg \
          - d_v_pos * dwdv_pos - d_v_neg * dwdv_neg)

    # flatten for output
    y_ = y_.flatten()

    return y_

# solve_ivp(fun, t_span, y0, method='RK45', t_eval=None, **)
ivp = wxv_ivp.flatten()
sol = solve_ivp(fokker_planck_system, [0, 2], ivp, dense_output=True)

# save the sol for later use
np.save('./sol.npy', sol)
