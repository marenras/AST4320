import numpy as np
import matplotlib.pyplot as plt


omega_m = 1.0
omega_l = 0.0


def solution(omega_m, omega_l):
    H_0 = 2.269e-18
    n = 200
    a0 = np.exp(-3); a1 = 1.0


    z1 = np.zeros(n+1)
    z2 = np.zeros(n+1)
    a, da = np.linspace(a0, a1, n+1, retstep = True)

    def H(a):
        return H_0*np.sqrt(omega_m*a**(-3) + omega_l)

    z1[0] = np.exp(-3.0)
    z2[0] = a[0]*H(a[0])


    def dzdt(a, z1, z2):
        return -2*H(a)*z2 + 3./2*H(a)**2*z1

    for i in range(n):
        dt = da/(a[i] * H(a[i]))
        z2[i+1] = z2[i] + dzdt(a[i], z1[i], z2[i]) * dt
        z1[i+1] = z1[i] + z2[i] * dt

    return np.log(a), np.log(z1)


a1, z1 = solution(1.0, 0.0)
a2, z2 = solution(0.3,0.7)
a3, z3 = solution(0.8,0.2)

plt.plot(a1, z1, 'r')
plt.plot(a2, z2, 'b')
plt.plot(a3, z3,'g')
plt.show()
