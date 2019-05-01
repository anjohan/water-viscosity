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

press = np.column_stack((pxy, pxz, pyz)) * 1e5

skip = 20000
interval_length = 20000
interval_distance = 2000
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
    interval_length,
    interval_distance,
    press.shape[0],
)
print(f"Loop time: {time.time() - t0}")

mean_viscosity *= V / (k_B * T) * 0.5e-15
stddev_viscosity *= V / (k_B * T) * 0.5e-15

t = np.arange(interval_length) * 5e-4
plt.style.use("ggplot")
plt.subplots(2, 1, sharex=True)
plt.subplot(2, 1, 1)
plt.xlabel("t [ps]")
plt.ylabel(r"$\langle P_{ij}(t)P_{ij}(0)\rangle$ [Pa$^2$]")
plt.plot(t, mean_autocorr)
plt.subplot(2, 1, 2)
plt.xlabel("t [ps]")
plt.ylabel(
    r"$\eta(t) = \frac{V}{k_\mathrm{B}T}\int_0^t\langle P_{ij}(t')P_{ij}(0)\rangle\,\mathrm{d}t'$ [Pa$\,$s]"
)
# plt.errorbar(t, mean_viscosity, stddev_viscosity)
plt.plot(t, mean_viscosity)
plt.tight_layout()
plt.savefig("data/corr_press.png", bbox_inches="tight")
plt.show()

# np.savetxt("data/corr.dat", np.column_stack((t, autocorr)))
# np.savetxt("data/viscosity.dat", np.column_stack((t, viscosity)))
