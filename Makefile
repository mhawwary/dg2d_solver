
# Specifying some environmental variables for Linux, note this can be done in the shell prompt

COMP	= GCC
TECIO	= NO
CODE	= DEBUG
OPENMP	= NO

# Specifing Standard Variables:
CXX	= g++ -std=gnu++11 #-pedantic-errors # c++ gcc compiler
CXXFLAGS=       # C++ compiler flags
LDLFLAGS=	# linker flags
CPPFLAGS=	# c/c++ preprocessor flags

OPTS	= 	# optimization flags and other options

# Includes

OPTS	+= -I include
#OPTS    += -I /home/mhawwary/Libraries/Eigen3.3.2/Eigen
#OPTS += -I /home/mhawwary/work/hpMusic/contrib/eigen/Eigen

ifeq ($(TECIO),YES)
	OPTS += -I $(TECIO_DIR)/tecsrc
endif

ifeq ($(CODE),RELEASE)
	ifeq ($(COMP),GCC)
		OPTS	+= -O3 
	endif
	
	ifeq ($(COMP),INTEL)
		OPTS	+= -xHOST -fast
	endif	
endif

ifeq ($(OPENMP),YES)	
	OPTS	+=  -lgomp -fopenmp 
endif

ifeq ($(CODE),DEBUG)
	ifeq ($(COMP),GCC)
		OPTS	+= -fbuiltin -g -Wall #-Werror
	endif
	
	ifeq ($(COMP),INTEL)
		OPTS	+= -xHOST -fast
	endif	
endif


# Source

SRC	= src/
OBJ	= obj/
BIN	= bin/
INC	= include/

vpath %.cpp src
vpath %.c src
vpath %.o   obj
vpath %.h src include  
vpath %.hpp src include 

# Objects
OBJS	= $(OBJ)DG2DFlow.o $(OBJ)SimData.o $(OBJ)SimCase.o $(OBJ)SimplyConnectGrid.o $(OBJ)Mesh.o $(OBJ)general_tools.o $(OBJ)solver_tools.o $(OBJ)dg_2d_diffus_solver.o $(OBJ)ExplicitTimeSolver.o # objects 
INCLS	= 

# Compile

.PHONY: default help clean


default: all
help:	
	@echo 'help'

all: DG2DFlow.exe

DG2DFlow.exe: $(OBJS)
	$(CXX) $(OPTS) -o $(BIN)$@ $+


$(OBJ)%.o : %.cpp 
	$(CXX) $(OPTS) -c -o $@ $<

$(OBJ)%.o : %.c 
	$(CXX) $(OPTS) -c -o $@ $<


$(OBJ)DG2DFlow.o:   DG2DFlow.cpp 
$(OBJ)SimData.o:   SimData.hpp SimData.cpp
$(OBJ)SimCase.o:   SimCase.hpp SimCase.cpp
$(OBJ)SimplyConnectGrid.o: SimplyConnectGrid.hpp SimplyConnectGrid.cpp
$(OBJ)Mesh.o: Mesh.hpp Mesh.cpp
$(OBJ)general_tools.o: general_tools.h general_tools.cpp
$(OBJ)solver_tools.o: solver_tools.h solver_tools.c
$(OBJ)dg_2d_diffus_solver.o: dg_2d_diffus_solver.hpp dg_2d_diffus_solver.cpp
$(OBJ)ExplicitTimeSolver.o: ExplicitTimeSolver.hpp ExplicitTimeSolver.cpp

clean:
	rm -f ./$(OBJ)*.o ./$(BIN)*.exe  
	@echo  removing all object and executable files

clean_temp:
	rm -f ./input/*.*~ ./$(OBJ)*.*~ ./$(BIN)*.*~ ./$(SRC)*.*~ ./$(INC)*.*~ ./*.*~ ./*~
	rm -f ./Results/*.*~ ./Results/*/*.*~ ./Results/*/*/*.*~ ./Results/*/*/*/*.*~ 
	

plot:
	python python_tools/DGsolplot_reader.py -f ./input/python_input.in 

plot_compare:
	python python_tools/DGplot_compare_Beta.py -f ./input/python_input_compare_Beta.in 

plot_errors:
	python python_tools/DGplot_errors.py -f ./input/python_input_errors.in 

plot_error_analys:
	python python_tools/DGplot_error_analysis.py -f ./input/python_input_error_analysis.in 


