solvers:
	rm -f solvers*.so
	f2py --f90flags="-Wall -Wextra" -m solvers -c solvers.f90
