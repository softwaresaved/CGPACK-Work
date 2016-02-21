#$Id: Makefile-archer 228 2016-02-19 14:23:28Z mexas $

#*robodoc*b* CGPACK/Makefile-archer
#  NAME
#    Makefile-archer
#  SYNOPSIS

FC=		ftn

#  PURPOSE
#    Build/install CGPACK with Cray compiler
#  AUTHOR
#    Anton Shterenlikht
#  COPYRIGHT
#    See CGPACK_Copyright
#  SOURCE

FFLAGS=		-c -eacFn -dm -m2 -rl
CGPACK=		cgpack
MYLIB=		lib$(CGPACK).a
LIBDIR=		$(HOME)/lib

# level 1

L1_SRC=		cgca_m1co.f90
L1_OBJ=		$(L1_SRC:.f90=.o)

# level 2 modules and submodules

L2_SRC_MOD=	cgca_m2alloc.f90 cgca_m2gb.f90 cgca_m2geom.f90 \
		cgca_m2glm.f90 cgca_m2hx.f90 cgca_m2lnklst.f90 \
		cgca_m2out.f90 cgca_m2pck.f90 cgca_m2phys.f90 \
		cgca_m2red.f90 cgca_m2rnd.f90 cgca_m2rot.f90 \
		cgca_m2stat.f90
L2_SRC_SUBMOD=	m2out_sm1.f90 m2out_sm2_mpi.f90
L2_SRC=		$(L2_SRC_MOD) $(L2_SRC_SUBMOD)
L2_OBJ=		$(L2_SRC:.f90=.o)

# level 3 modules and submodules 

L3_SRC_MOD=	cgca_m3clvg.f90 cgca_m3gbf.f90 cgca_m3nucl.f90 \
		cgca_m3pfem.f90 cgca_m3sld.f90
L3_SRC_SUBMOD=	m3clvg_sm1.f90 m3clvg_sm2.f90 m3clvg_sm3.f90 \
		m3sld_sm1.f90
L3_SRC=		$(L3_SRC_MOD) $(L3_SRC_SUBMOD)
L3_OBJ=		$(L3_SRC:.f90=.o)

# level 4

L4_SRC_MOD=	cgca_m4fr.f90
L4_OBJ=		$(L4_SRC_MOD:.f90=.o)

# top level

LTOP_SRC_MOD=	cgca.F90
LTOP_OBJ=	$(LTOP_SRC_MOD:.F90=.o)

ALL_OBJ=	$(L1_OBJ) $(L2_OBJ) $(L3_OBJ) $(L4_OBJ) $(LTOP_OBJ)

.SUFFIXES: .F90 .f90 .o

all: $(ALL_OBJ) $(MYLIB)

.f90.o:
	$(FC) $(FFLAGS) $<

.F90.o:
	$(FC) $(FFLAGS) $<

# object dependencies

$(L2_OBJ): $(L1_OBJ)
$(L3_OBJ): $(L2_OBJ)
$(L4_OBJ): $(L3_OBJ)
$(LTOP_OBJ): $(L4_OBJ)

# Submodule dependencies. Make each submodule depend only on its own
# parent module. So have to use exact names here.
# level 2
m2out_sm1.o: cgca_m2out.o

# level 3
m3clvg_sm1.o m3clvg_sm2.o m3clvg_sm3.o: cgca_m3clvg.o
m3sld_sm1.o: cgca_m3sld.o

$(MYLIB): $(ALL_OBJ)
	@if [ -e $(MYLIB) ]; then \
		rm $(MYLIB); \
	fi
	ar -r $(MYLIB) $(ALL_OBJ)

install: $(MYLIB)
	@if [ -e $(LIBDIR)/$(MYLIB) ]; then \
		echo $(LIBDIR)/$(MYLIB) already exists; \
		echo run \"make deinstall\" first; \
		exit 1; \
	fi
	cp $(MYLIB) $(LIBDIR)

deinstall:
	cd $(LIBDIR) && rm $(MYLIB)

clean:
	rm $(ALL_OBJ) $(MYLIB) *.lst

#*roboend*
