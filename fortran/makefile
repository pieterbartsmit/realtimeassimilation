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
opt  = -Og -g -fcheck=all -fbacktrace -fopenmp
file = modPar.F90 modIO.F90 modTools.f90 modRayIntegration.f90 raytracing.f90

all:
	$(CC) $(flag) -shared  -fPIC -fbacktrace $(opt) $(file) -o raytracing.lib

clean:
	cd ./compile/
