FC= ifort
OBJS=  param.f90  sub.f90  main.f90

build:  
	$(FC) -o grav $(OBJS) 

#.SUFFIXS : .f90

clean:
	rm -f *.o *.mod fastslow
run:
	make build
	./grav
