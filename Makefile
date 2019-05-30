all: data/corr_press.png

data/corr_press.png: analysis.py data/log.viscosity autocorr.cpython-37m-x86_64-linux-gnu.so
	python $<

data/log.viscosity: in.viscosity data/data.minimised
	lmp -k on g 1 -sf kk -pk kokkos newton on neigh half binsize 7.5 -in $<

data/data.minimised: in.minimise data/data.water
	mpirun lmp -in $<

data/data.water: run_packmol.py xyz.water
	python $<

%.cpython-37m-x86_64-linux-gnu.so: %.f90
	f2py -c -m $* -DF2PY_REPORT_ON_ARRAY_COPY --opt="-fcheck=all -O3 -fopenmp" $<
