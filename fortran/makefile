#
#-----------------------------------------------------------
#              Makefile for the SO Wave Theory Files
#-----------------------------------------------------------
#
# == Variables ==
#
# Compiler 
CC = gfortran
#
# Workdirectory
fdir = ./
wdir = ./

#
# Flags
#
flag = -I $(wdir) -J $(wdir)
opt  = -o3 -ffast-math  -fopenmp -march=core2
file = modPar.F90 modIO.f90 modTools.f90 modRayIntegration.f90 RayTracing.f90

all:
	$(CC) $(flag) -shared  -static $(opt) $(file) -o raytracing.lib

clean:
	cd ./compile/
