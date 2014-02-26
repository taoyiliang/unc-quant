## Put all macros needed below in Makefile.include, except for those
## defined in the nox export files.  Note some of these may be defined
## in environment variables.


#
# include file from NOX.
#
TRILINOS_DIR=$(HOME)/packages/trilinos-11.4.1-install-gcc4.6.3
include $(TRILINOS_DIR)/include/Makefile.export.NOX

  TRILINOS_INSTALL_DIR=$(TRILINOS_DIR)
  COMPILE_FLAGS = $(NOX_CXX_FLAGS) $(NOX_TPL_INCLUDES) \
	  	-I$(TRILINOS_INSTALL_DIR)/include -I./include 
  CXX=$(NOX_CXX_COMPILER) #need to redefine these so that the implicit rule .o.C works correctly
  FC=$(NOX_FORTRAN_COMPILER) #need to redefine these so that the implicit rule .o.C works correctly
  LD=$(NOX_CXX_COMPILER)
  CPPFLAGS=$(COMPILE_FLAGS) #-DMPICH_IGNORE_CXX_SEEK
  MY_CXX_FLAGS=$(NOX_CXX_FLAGS) $(PBFLAG)

##
## Build linker flags
##
LINK_FLAGS = -L$(TRILINOS_INSTALL_DIR)/lib $(NOX_LIBRARIES) $(NOX_TPL_LIBRARIES) $(NOX_EXTRA_LD_FLAGS) -lgfortran 

##
## Build application code
##
SRC 	:= $(wildcard src/*.C) $(wildcard src/*/*.C) $(wildcard src/*/*/*.C)
FSRC	:= $(wildcard src/fortran/*.f)

objects	:= $(patsubst %.C, %.o, $(SRC))
fobjects:= $(patsubst %.f, %.o, $(FSRC))

include_DIR := $(NOX_TPL_INCLUDES) -I$(TRILINOS_INSTALL_DIR)/include -I./include 


TwoDProblem: $(objects) $(fobjects) .depend
	$(LD) $(MY_CXX_FLAGS) -o TwoDProblem $(objects) $(fobjects) $(LINK_FLAGS)

%.$(obj-suffix) : %.C
	@echo "Compiling C "$<"..."
	$(CXX) $(COMPILE_FLAGS)  -c $< -o $@

%.$(fobj-suffix) : %.f
	@echo "Compiling fortran "$<"..."
	$(FC) $(COMPILE_FLAGS)  -c $< -o $@
clean:
	rm -f $(objects)  $(fobjects) *.o TwoDProblem  *~
	rm .depend
info:
	@echo "TRILINOS_INSTALL_DIR " $(TRILINOS_INSTALL_DIR)
	@echo "COMPILE_FLAGS " $(COMPILE_FLAGS)
	@echo "LINK_FLAGS " $(LINK_FLAGS)
	@echo "SRC: " $(SRC)
	@echo "OBJ: " $(OBJS)
	@echo "CXX: " $(CXX)
	@echo "LD: " $(LD)
	@echo "NOX_CXX_FLAGS: " $(NOX_CXX_FLAGS)
	@echo "NOX_LIBRARIES: " $(NOX_LIBRARIES)
	@echo "NOX_TPL_LIB: " $(NOX_TPL_LIBRARIES)
	@echo "NOX_EXTRA_LD_FLAGS: " $(NOX_EXTRA_LD_FLAGS)
	@echo "debug_objs: " $(debug_objs)
	@echo "MY_CXX_FLAG: " $(MY_CXX_FLAG)

.PHONY: method

# include the dependency list
include .depend


#
# Dependencies
#
.depend:
	@perl make_dependencies.pl -I. $(include_DIR) "-S\$$(obj-suffix)" $(SRC) > .depend
	@echo "Updated .depend"
