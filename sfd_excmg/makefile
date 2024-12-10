SRC_DIR = src
OBJ_DIR = obj

OBJECTS = obj/global_para.o obj/model_selec.o obj/file_opt.o
OBJECTS+= obj/zagmg_mumps.o obj/zagmg.o
OBJECTS+= obj/div_correc.o obj/excmg.o
OBJECTS+= obj/solvers.o obj/fdfe_main.o
OBJECTS+= obj/sfd_excmg.o

BIN = sfd3d

F90 = ifort
MyCC = $(MPI_HOME)/intel64/bin/mpicc
C_COMPILE_FLAGS = -O3 -g -Wall -DHYPRE_TIMING
CFLAG = ${C_COMPILE_FLAGS} #${HYPREINCS}
FFLAG = -O3 -traceback -qopenmp
FCCFLAG = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64  -liomp5 -lpthread
mkllib=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64
mklinc=/opt/intel/oneapi/mkl/2022.0.2/include
HEDPS = $(OBJECTS:.o=.d)
MATHCLIB = -lm  -lirc -lsvml -lstdc++
MPILIBFLAG = -L$(MPI_HOME)/intel64/lib
#HYPRELIBFLAG = -L$(HYPRE_HOME)/lib -lHYPRE
ALLLIBS = $(MPILIBFLAG) #$(MATHCLIB)


export F90
export FFLAG
export CFLAG
vpath %.f90 $(SRC_DIR)
vpath %.f $(SRC_DIR)
vpath %.c $(SRC_DIR)
#------
$(BIN):$(OBJECTS)
	$(F90) -o $@ $^ -I$(mklinc) -L$(mkllib)  $(FCCFLAG)

#$(OBJ_DIR)/%.o: %.f90
#	$(F90) $(FFLAG) -c -o $@ $<
#$(OBJ_DIR)/%.o: %.f
#	$(F90) $(FFLAG) -c -o $@ $< 
#$(OBJ_DIR)/%.o: %.c
#	$(MyCC) $(CFLAG) -c -o $@ $<
#
#the dependencies of the modules

obj/global_para.o:src/global_para.f90
	$(F90) $(FFLAG) -c -o $@ $<
obj/model_selec.o:src/model_selec.f90 src/global_para.f90
	$(F90) $(FFLAG) -c -o $@ $< 
obj/file_opt.o:src/file_opt.f90 src/global_para.f90
	$(F90) $(FFLAG) -c -o $@ $<
obj/zagmg_mumps.o:src/zagmg_mumps.f90
	$(F90) $(FFLAG) -c -o $@ $<
obj/zagmg.o:src/zagmg.f90
	$(F90) $(FFLAG) -c -o $@ $<
obj/div_correc.o:src/div_correc.f90 src/global_para.f90 src/model_selec.f90
	$(F90) $(FFLAG) -c -o $@ $<
obj/excmg.o:src/excmg.f90 src/global_para.f90
	$(F90) $(FFLAG) -c -o $@ $<
obj/solvers.o:src/solvers.f90 src/div_correc.f90 src/global_para.f90 
	$(F90) $(FFLAG) -c -o $@ $<
obj/fdfe_main.o:src/fdfe_main.f90 src/global_para.f90 src/solvers.f90 src/model_selec.f90 src/div_correc.f90\
 src/file_opt.f90 src/excmg.f90 src/zagmg.f90 src/zagmg_mumps.f90
	$(F90) $(FFLAG) -c -o $@ $<
obj/sfd_excmg.o:src/sfd_excmg.f90 src/fdfe_main.f90 src/global_para.f90 src/model_selec.f90 src/file_opt.f90
	$(F90) $(FFLAG) -c -o $@ $<

-include $(HDEPS)


clean:
	rm -f obj/*.o
	rm -f obj/*.d
	rm -f *.mod
            
