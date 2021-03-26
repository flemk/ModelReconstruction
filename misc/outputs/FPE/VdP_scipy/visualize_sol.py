import numpy as np
import matplotlib.pyplot as plt

bins = 50

print('loading start')
sol = np.load('./sol.npy', allow_pickle=True)
print('loading end')
print('---')

t = np.linspace(0, 5, 300)
z = sol.any(0).sol(t)

z = z.T.reshape((300, bins, bins))

for t in [0, 25, 35, 50, 100, 150, 160, 170, 180, 200, 250, 299]:
    fig = plt.figure()
    o = plt.contourf(z[t])
    fig.colorbar(o)
    plt.title('t = %s' % (str(t)))
    plt.savefig('%s.noupload.png' % (str(t)))
    plt.clf()
