import numpy as np
import matplotlib.pyplot as plt

R = 1.0
i = np.sqrt(-1+0j)

def W(x):
    """
    Top-hat smoothing function.
    """
    return (abs(x) < R)*1.0

def W_tilde(k):
    """
    Fourier transformation of W.
    """
    #return -1.0/(i*k) * np.exp(-i*k*R) + 1.0/(i*k) * np.exp(i*k*R)
    return 2.0/k * np.sin(k*R)


def FWHM(x_values, f_values):
    """
    Calculates the full width at half maximum.
    """
    half_max = 0.5*np.max(f_values)
    x_fwhm = 0
    for i in range(len(x_values)):
        if f_values[i] < half_max:
            x_fwhm = x_values[i]
        else:
            break
    return (abs(x_fwhm)*2)


N = 1000
x = np.linspace(-2,2, N)
k = np.linspace(-20,20, N)

print("Full width at half maximum: %.2f" %(FWHM(k,W_tilde(k))))

plt.figure()
plt.plot(x, W(x), color='magenta')
plt.title(r'Top-hat smoothing function $W$')
plt.xlabel('x')
plt.grid()

plt.figure()
plt.plot(k, W_tilde(k), color='magenta')
plt.title(r'Fourier conjugate $\tilde{W}$ of top-hat smoothing function')
plt.xlabel('k')
plt.grid()
plt.show()
