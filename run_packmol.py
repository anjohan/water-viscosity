import os
import sys
import numpy as np

density = 0.0335
N = 2000
L = (N / density) ** (1.0 / 3)

packmol_script = f"""
tolerance 2.0

filetype xyz
output data/xyz.water

nloop 10000

structure xyz.water
    number {N}
    inside box 0. 0. 0. {L} {L} {L}
end structure
"""

with open("data/packmol_script.inp", "w") as outfile:
    outfile.write(packmol_script)

packmol_return_code = os.system("packmol < data/packmol_script.inp")
if packmol_return_code != 0:
    sys.exit(packmol_return_code)

data = np.loadtxt("./data/xyz.water", unpack=True, skiprows=2)
number_of_atoms = len(data[0])

new_data = np.column_stack((np.arange(1, number_of_atoms + 1), *data))

np.savetxt(
    "./data/data.water",
    new_data,
    fmt="%d %d %g %g %g",
    header=f"""created by me

{3*N} atoms
3 atom types

0 {L} xlo xhi
0 {L} ylo yhi
0 {L} zlo zhi

Masses

1 28.08
2 15.9994
3 1.00794

Atoms
""",
    comments="",
)
