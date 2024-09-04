import numpy as np
import matplotlib.pyplot as plt

def autocorrelation(v, dt):
    n = len(v)
    result = np.correlate(v, v, mode='full')
    result = result[result.size // 2:] / np.arange(n, 0, -1)
    return result * dt  

dt = 0.001  

data_mag = np.loadtxt('velMAG.dat')
Vx_mag = data_mag[:, 0]
Vy_mag = data_mag[:, 1]
V_mag = np.sqrt(Vx_mag**2 + Vy_mag**2)

data_nm = np.loadtxt('velNM.dat')
Vx_nm = data_nm[:, 0]
Vy_nm = data_nm[:, 1]
V_nm = np.sqrt(Vx_nm**2 + Vy_nm**2)

Av_mag = autocorrelation(V_mag, dt)
Av_nm = autocorrelation(V_nm, dt)

tau = np.arange(len(Av_mag)) * dt

plt.plot(tau, Av_mag, label='AV (Con Magnus)')
plt.plot(tau, Av_nm, label='AV (Sin Magnus)')
plt.ylim(0, 0.2)
plt.xlabel(r'$\tau$ (s)')
plt.ylabel(r'$A_v(\tau)$')
plt.legend()
plt.title('Función de Autocorrelación de la Velocidad')
plt.savefig("Corr.png")

#prueba

