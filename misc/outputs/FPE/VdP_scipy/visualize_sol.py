import numpy as np
import matplotlib.pyplot as plt

bins = 200

print('loading start')
sol = np.load('./sol.npy', allow_pickle=True)
print('loading end')
print('---')

t = np.linspace(0, 5, 300)
z = sol.any(0).sol(t)

z = z.T.reshape((300, bins, bins))


plt.contourf(z[0])
plt.title('t = 0')
plt.savefig('t0.png')
plt.clf()

plt.contourf(z[50])
plt.title('t = 50')
plt.savefig('t50.png')
plt.clf()

plt.contourf(z[100])
plt.title('t = 100')
plt.savefig('t100.png')
plt.clf()

plt.contourf(z[150])
plt.title('t = 150')
plt.savefig('t150.png')
plt.clf()

plt.contourf(z[200])
plt.title('t = 200')
plt.savefig('t200.png')
plt.clf()

plt.contourf(z[250])
plt.title('t = 250')
plt.savefig('t250.png')
plt.clf()

plt.contourf(z[-1])
plt.title('t = -1')
plt.savefig('t-1.png')
plt.clf()
