FC= ifort
OBJS=  param.f90 solver.f90  sub.f90  main.f90

build:  
	$(FC) -o chase  $(OBJS) 

#.SUFFIXS : .f90

clean:
	rm -f *.o *.mod fastslow
run:
	make build
	./chase
