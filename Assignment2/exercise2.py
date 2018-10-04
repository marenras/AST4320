import numpy as np
import matplotlib.pyplot as plt


def random_walk():
    """
    Function for doing random walk by reducing S_c with a
    small number epsilon for each step.
    """
    # Creating initial values
    var = 0.5*1e-4
    S_c = [(np.pi/var)**(1.0/4)]
    delta = [np.random.normal(scale=np.sqrt(var))]
    eps = 0.1

    # Doing the random walk untill S_c reaches 1.
    i = 1
    while S_c[-1] >= 1.0:
        S_c.append(S_c[i-1] - eps)
        var_new = np.pi/S_c[i]**4
        beta = np.random.normal(scale=(np.sqrt(var_new-var)))
        delta.append(delta[i-1] + beta)
        var = var_new
        i+=1

    return delta, S_c


def gauss_dist(delta, S_c):
    """
    Gaussian distribution for overdensity smoothed with
    scale S_c, corresponding to mass M, for variance
    sigma^2(M).
    """
    var = (np.pi/(S_c**4))
    return 1.0/(np.sqrt(2*np.pi*var))*np.exp(-delta**2/(2*var))


def gauss_dist_nc(delta, S_c):
    """
    Gaussian distribution for overdensity smoothed with
    scale S_c, corresponding to mass M, for variance
    sigma^2(M), when the chains never cross delta_crit=1.0.
    """
    var = (np.pi/(S_c**4))
    delta_crit = 1.0
    return 1.0/(np.sqrt(2*np.pi*var))*(np.exp(-delta**2/(2*var)) - \
    np.exp(- (2*delta_crit - delta)**2/(2*var)))



deltas = []
S_cs = []
deltas_nc = []
S_cs_nc = []

# Doing the random walk 10^5 times and saving the last values
N = int(1e5)
for j in range(N):
    delta, S_c = random_walk()
    deltas.append(delta[-1])
    S_cs.append(S_c[-1])
    if delta[-1] < 1:
        deltas_nc.append(delta[-1])
        S_cs_nc.append(S_c[-1])


# Creating array with delta values between the max and
# min delta from the random walk.
delta_array = np.linspace(np.min(deltas),np.max(deltas),len(deltas))
delta_nc_array = np.linspace(np.min(deltas_nc),np.max(deltas_nc),len(deltas_nc))


plt.figure()
plt.title(r'Random walk for $\delta$-values giving $S_c = 1.0$')
plt.hist(deltas, bins=40, density = True, color='teal')
plt.plot(delta_array, gauss_dist(delta_array, np.array(S_cs)), color='magenta')

plt.figure()
plt.title(r'Random walk for $\delta$-values never crossing $\delta_{crit} = 1.0$')
plt.hist(deltas_nc, bins=40, density = True, color='teal')
plt.plot(delta_nc_array, gauss_dist_nc(delta_nc_array, np.array(S_cs_nc)), color='magenta')
plt.show()
