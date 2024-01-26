run:quadrature.o main.o
	gfortran -o $@ $^ -fopenmp

%.o:%.f90
	gfortran -c -O2 -o $@ $< -fopenmp

clean:
	@rm -f run *.o *.mod *.vtu
