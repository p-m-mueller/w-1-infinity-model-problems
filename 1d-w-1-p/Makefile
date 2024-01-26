all: debug run

debug: main.f90
	gfortran -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero -O0 -Wall --check=all -pedantic $< -llapack

release: main.f90 
	gfortran -O2 $< -llapack

run: a.out
	./a.out

clean:
	@rm -v a.out *.o *.dat run

plot:
	@gnuplot plot-results.gnuplot

animate:
	@gnuplot animate-results.gnuplot
#	@ffmpeg -y -framerate 10 -pattern_type glob -i 'out/animation/*.png' -c:v ffv1 out/animation.avi
	@ffmpeg -y -framerate 10 -pattern_type glob -i 'out/animation/*.png' -c:v libx265 out/animation.mp4