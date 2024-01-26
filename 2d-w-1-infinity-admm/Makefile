run:quadrature.o main.o
	gfortran -o $@ $^ -fopenmp

%.o:%.f90
	gfortran -c -O3 -o $@ $< -fopenmp

clean:
	rm -vf run *.o *.mod *.vtu *.dat
