import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import time

dt = 0.01
T = 100000 * np.pi
n = int(T / dt)
sqrtdt = np.sqrt(dt)

def f(t, x):
    '''Van-der-Pol oscillator
    '''
    epsilon = 2
    g = 3
    y = [0, 0]
    # 1. deterministic
    y[0] = x[1]
    y[1] = (epsilon - x[0] ** 2) * x[1] - x[0]
    # 2. stochastic
    # noise only in y[1]
    y[1] += sqrtdt * np.random.randn() * g
    return y

# get runtime
start_time = time.time()

# solve now
sol = solve_ivp(f, [0, T], [0, 1], dense_output=True)
t = np.linspace(0, T, n)
x, v = sol.sol(t)

# runtime done
runtime = time.time() - start_time
print('time elapsed: %s' %(str(runtime)))

# save those
np.save('VdP_series_xv.npy', (x, v))
np.save('VdP_scipy_sol.npy', sol)

# save a phaseplot
plt.scatter(x, v, s=1, label='phaseplot')
plt.title('Phaseplot of VdP white noised')
plt.xlabel('x')
plt.ylabel('v')
plt.savefig('VdP_phaseplot.png')
