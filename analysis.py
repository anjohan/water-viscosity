import matplotlib.pyplot as plt
import numpy as np
from tqdm import trange
from logplotter import find_data
import time
import autocorr

calc_autocorr = autocorr.mod_autocorr.calc_autocorr

density = 0.0335
N = 2000
L = (N / density) ** (1.0 / 3) * 1e-10
V = L ** 3
T = 300.0
k_B = 1.38e-23

data = find_data("data/log.viscosity")

step = data["Step"]
pxy = np.asarray(data["Pxy"])
pxz = np.asarray(data["Pxz"])
pyz = np.asarray(data["Pyz"])

print(pxy.mean(), pxz.mean(), pyz.mean())

press = np.column_stack((pxy, pxz, pyz)) * 1e5

skip = 20000
interval_length = 20000
interval_distance = 200
num_intervals = 3 * (press.shape[0] - interval_length + 1) // interval_distance

press = np.asfortranarray(press[skip:, :])

mean_autocorr, mean_viscosity, stddev_viscosity = np.zeros((3, interval_length))
autocorrs = np.zeros((interval_length, num_intervals), order="F")
viscosities = np.zeros((interval_length, num_intervals), order="F")

t0 = time.time()
calc_autocorr(
    press,
    autocorrs,
    viscosities,
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
viscosities *= V / (k_B * T) * 0.5e-15

t = np.arange(interval_length) * 5e-4
plt.figure(figsize=(6, 9))
plt.style.use("ggplot")
ax1 = plt.subplot(3, 1, 1)
plt.ylabel(r"$\langle P_{ij}(t)P_{ij}(0)\rangle$ [Pa$^2$]")
for i in range(0, autocorrs.shape[1], autocorrs.shape[1] // 100):
    plt.plot(t, autocorrs[:, i])
plt.plot(t, mean_autocorr, "k-", lw=4.0)
plt.subplot(3, 1, 2, sharex=ax1)
# plt.xlabel("t [ps]")
plt.ylabel(
    r"$\eta(t) = \frac{V}{k_\mathrm{B}T}\int_0^t\langle P_{ij}(t')P_{ij}(0)\rangle\,\mathrm{d}t'$ [Pa$\,$s]"
)
# plt.errorbar(t, mean_viscosity, stddev_viscosity)
for i in range(0, viscosities.shape[1], viscosities.shape[1] // 100):
    plt.plot(t, viscosities[:, i])
plt.plot(t, mean_viscosity, "k-", lw=4.0)
plt.subplot(3, 1, 3, sharex=ax1)
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

np.savetxt("data/corr.dat", np.column_stack((t, mean_autocorr)))
np.savetxt("data/viscosity.dat", np.column_stack((t, mean_viscosity, stddev_viscosity)))

estimate_step = 2000
print(mean_viscosity[estimate_step], stddev_viscosity[estimate_step])
