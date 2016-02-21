#$Id: Makefile-archer 149 2015-10-27 15:10:14Z mexas $

#*robodoc*b* tests/Makefile-tests-archer
#  NAME
#    Makefile-tests-archer
#  SYNOPSIS

FC=		ftn

#  PURPOSE
#    Build CGPACK tests on HECToR with Cray compilers.
#  AUTHOR
#    Anton Shterenlikht
#  COPYRIGHT
#    See CGPACK_Copyright
#  SOURCE

#FFLAGS=	-eacn -Rb -I. -I$(LIBDIR) 
#FFLAGS=	-eacn -Rb -m1 -rl -I. -I$(LIBDIR) 
FFLAGS=         -c -eacFn -dm -m1 -rl -I. -I$(LIBDIR)

CGLIB=		cgpack
LIBDIR=		$(HOME)/lib
LIBNAME=	lib$(CGLIB)

MODSRC=		testaux.f90
MODOBJ=		$(MODSRC:.f90=.o)
MODLST=		$(MODSRC:.f90=.lst)
MOD_CLEAN=	$(MODLST) $(MODOBJ)

SRC= \
testAAA.f90 testAAB.f90 testAAC.f90 testAAD.f90 testAAE.f90 \
testAAF.f90 testAAG.f90 testAAH.f90 testAAI.f90 testAAJ.f90 \
testAAK.f90 testAAL.f90 testAAM.f90 testAAN.f90 testAAO.f90 \
testAAP.f90 testAAQ.f90 testAAR.f90 testAAS.f90 testAAT.f90 \
testAAU.f90 testAAV.f90 testAAW.f90 testAAX.f90 testAAY.f90 \
testAAZ.f90 \
testABA.f90 testABB.f90 testABC.f90 testABD.f90 testABE.f90 \
testABF.f90 testABG.f90 testABH.f90 testABI.f90 testABJ.f90 \
testABK.f90 testABL.f90 testABM.f90 testABN.f90 testABO.f90 \
testABP.f90 testABQ.f90 testABR.f90 testABS.f90 testABT.f90 \
testABU.f90 testABV.f90 testABW.f90 testABX.f90 testABY.f90 \
testABZ.f90

NON_COARRAY_SRC=test_gc.f90
NON_COARRAY_EXE=$(NON_COARRAY_SRC:.f90=.xnonca)

OBJ=		${SRC:.f90=.o}
EXE=		${SRC:.f90=.x} ${NON_COARRAY_EXE}
LST=		$(SRC:.f90=.lst)

CAMPI_SRC=	ca_mpi.f90
CAMPI_EXE=	$(CAMPI_SRC:.f90=.x)
CAMPI_LST=	$(CAMPI_SRC:.f90=.lst)
CAMPI_CLEAN=	$(CAMPI_EXE) $(CAMPI_LST)

ALL_CLEAN=	$(MOD_CLEAN) $(OBJ) $(EXE) $(LST) $(CAMPI_CLEAN)

.SUFFIXES: .f90 .o .x .mod .xnonca

LOGIN_NODE_FLAGS=-h cpu=x86-64 -eacn -Rb
LIB=-L$(LIBDIR) -l$(CGLIB)

all: $(OBJ) $(EXE) $(CAMPI_EXE)

.f90.o:
	$(FC) -c $< $(FFLAGS)

.o.x:
	$(FC) $< -o $@ $(MODOBJ) $(LIB)

.f90.xnonca:
	$(FC) -o $@ $< $(LOGIN_NODE_FLAGS)

$(OBJ): $(MODOBJ) $(LIBDIR)/$(LIBNAME).a
$(MODOBJ): $(LIBDIR)/$(LIBNAME).a
 
clean:
	\rm $(ALL_CLEAN)

#*roboend*
