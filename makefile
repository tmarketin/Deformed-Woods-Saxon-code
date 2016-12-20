FC = gfortran

OBJ = params.o func.o eigen.o phys.o wf.o dws.o

%.o: %.f90
		$(FC) -c -o $@ $<

defomedWS: $(OBJ)
	$(FC) -o run_deformedWS $(OBJ)

