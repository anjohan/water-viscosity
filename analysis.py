import matplotlib.pyplot as plt
import numpy as np
from tqdm import trange
from logplotter import find_data
import time
import autocorr

calc_autocorr = autocorr.mod_autocorr.calc_autocorr

density = 0.0335
N = 6000
L = (N / density) ** (1.0 / 3) * 1e-10
V = L ** 3
T = 300.0
k_B = 1.38e-23

t0 = time.time()
data = find_data("15ns_data/log.viscosity")
print(f"Read time: {time.time() - t0}")

step = data["Step"]
pxy = np.asarray(data["Pxy"])
pxz = np.asarray(data["Pxz"])
pyz = np.asarray(data["Pyz"])

press = np.column_stack((pxy, pxz, pyz)) * 1e5

skip = 20000
interval_length = 4000
interval_distance = 600
num_intervals = 3 * (press.shape[0] - interval_length + 1) // interval_distance

press = np.asfortranarray(press[skip:, :])

mean_autocorr, mean_viscosity, stddev_viscosity = np.zeros((3, interval_length))

t0 = time.time()
calc_autocorr(
    press,
    mean_autocorr,
    mean_viscosity,
    stddev_viscosity,
    window_length=interval_length,
    window_distance=interval_distance,
    num_samples=press.shape[0],
)
print(f"Loop time: {time.time() - t0}")

mean_viscosity *= V / (k_B * T) * 0.5e-15
stddev_viscosity *= V / (k_B * T) * 0.5e-15 / np.sqrt(num_intervals)

t = np.arange(interval_length) * 5e-4
plt.figure(figsize=(6, 9))
plt.style.use("ggplot")
ax1 = plt.subplot(2, 1, 1)
plt.ylabel(r"$\langle P_{ij}(t)P_{ij}(0)\rangle$ [Pa$^2$]")
plt.plot(t, mean_autocorr)
plt.subplot(2, 1, 2, sharex=ax1)
plt.xlabel("t [ps]")
plt.ylabel(
    r"$\eta(t) = \frac{V}{k_\mathrm{B}T}\int_0^t\langle P_{ij}(t')P_{ij}(0)\rangle\,\mathrm{d}t'$ [Pa$\,$s]"
)
plt.plot(t, mean_viscosity)
plt.fill_between(
    t, mean_viscosity - stddev_viscosity, mean_viscosity + stddev_viscosity, alpha=0.2
)
plt.tight_layout()
plt.savefig("data/corr_press.png", bbox_inches="tight")
plt.show()

np.savetxt(
    "/home/anders/master/data/viscosity/corr_367.dat",
    np.column_stack((t, mean_autocorr)),
)
np.savetxt(
    "/home/anders/master/data/viscosity/viscosity_367.dat",
    np.column_stack((t, 1e3 * mean_viscosity, 1e3 * stddev_viscosity)),
)

estimate_step = 2000
print(mean_viscosity[estimate_step], stddev_viscosity[estimate_step])
