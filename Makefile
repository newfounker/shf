srcdat = dsyev.f
#		 calc_distance.f95\
#		 calc_angle.f95\
#		 calc_oop.f95\
#		 calc_dihedral.f95\
#		 calc_center.f95\
#		 calc_moments.f95\
#		 calc_rot.f95

outdat = dsyev.o
#		 calc_distance.o\
#		 calc_angle.o\
#		 calc_oop.o\
#		 calc_dihedral.o\
#		 calc_center.o\
#		 calc_moments.o\
#		 calc_rot.o

compile:
#	gfortran -std=f2003 -c $(srcdat)
	gfortran -g -llapack -lblas -std=f2003 -o hf-scf main.f95 #$(outdat)
