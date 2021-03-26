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
    if (t % 0.1 == 0):
        print('currently at t:%f' % (t))

    # calculate wxv "field" from 1d flattend array y
    wxv = y
    wxv = wxv.reshape((bins, bins))

    # apply neumann boundary conditions
    # copy paste from Andreas
    wxv[0, :] = (4.0 * wxv[1, :] - wxv[2, :]) / 3.
    wxv[-1, :] = (4.0 * wxv[-2, :] - wxv[-3, :]) / 3.
    wxv[:, -1] = (4.0 * wxv[:, -2] - wxv[:, -3]) / 3.
    wxv[:, 0] = (4.0 * wxv[:, 1] - wxv[:, 2]) / 3.
    # at the corners
    wxv[0, 0] = 4.0 * wxv[1, 1] - wxv[2, 2]
    wxv[-1, 0] = 4.0 * wxv[-2, 1] - wxv[-3, 2]
    wxv[0, -1] = 4.0 * wxv[1, -2] - wxv[2, -3]
    wxv[-1, -1] = 4.0 * wxv[-2, -2] - wxv[-3, -3]

    # flatten the coefficients
    d_x_pos_ = d_x_pos.flatten()
    d_x_neg_ = d_x_neg.flatten()
    d_v_pos_ = d_v_pos.flatten()
    d_v_neg_ = d_v_neg.flatten()

    # flatten the axis
    x_1_ = x_1.flatten()

    # receive first order upwind derivates
    dwdx_pos, dwdx_neg, dwdv_pos, dwdv_neg = first_order_upwind(wxv, dx)
    # set them to zero at the bounds
    '''dwdx_pos[0,:] = 0
    dwdx_pos[-1,:] = 0
    dwdx_pos[:,0] = 0
    dwdx_pos[:,-1] = 0

    dwdx_neg[0,:] = 0
    dwdx_neg[-1,:] = 0
    dwdx_neg[:,0] = 0
    dwdx_neg[:,-1] = 0

    dwdv_pos[0,:] = 0
    dwdv_pos[-1,:] = 0
    dwdv_pos[:,0] = 0
    dwdv_pos[:,-1] = 0

    dwdv_neg[0,:] = 0
    dwdv_neg[-1,:] = 0
    dwdv_neg[:,0] = 0
    dwdv_neg[:,-1] = 0'''
    # flatten those derivates
    dwdx_pos = dwdx_pos.flatten()
    dwdx_neg = dwdx_neg.flatten()
    dwdv_pos = dwdv_pos.flatten()
    dwdv_neg = dwdv_neg.flatten()

    # the second derivate of diffusion coefficient is constantly
    d2wdv2 = (np.roll(wxv, shift=(0, -1), axis=(1, 0)) - 2 * wxv + np.roll(wxv, shift=(0, 1), axis=(1, 0))) / (dx ** 2)
    # and flatten it aswell
    d2wdv2 = d2wdv2.flatten()

    # now apply for each point in grid y
    y_ = np.empty(len(y)) # initilize result array
    for i in range(len(y)):
        y_[i] = + (epsilon - x_1_[i] ** 2) * y[i] + (g / 2) * d2wdv2[i] \
                - d_x_pos_[i] * dwdx_pos[i] - d_x_neg_[i] * dwdx_neg[i] \
                - d_v_pos_[i] * dwdv_pos[i] - d_v_neg_[i] * dwdv_neg[i]

    return y_

# solve_ivp(fun, t_span, y0, method='RK45', t_eval=None, **)
ivp = wxv_ivp.flatten()
sol = solve_ivp(fokker_planck_system, [0, 2], ivp, dense_output=True)

# save the sol for later use
np.save('./sol.npy', sol)
