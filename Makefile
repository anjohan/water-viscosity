all: data/corr_press.png

data/corr_press.png: analysis.py data/log.viscosity autocorr.cpython-37m-x86_64-linux-gnu.so
	python $<

data/log.viscosity: in.viscosity data/data.water
	mpirun lmp -in $<

data/data.water: run_packmol.py
	python $<

%.cpython-37m-x86_64-linux-gnu.so: %.f90
	f2py -c -m $* -DF2PY_REPORT_ON_ARRAY_COPY --opt="-O3 -fopenmp" $<
