#
# Makefile for the DiscreteState project
#
# HOW TO USE:
# 1. Load the Intel environment first:
#    $ module load oneapi/compiler-intel-llvm/latest
# 2. To build the release version:
#    $ make
# 3. To build the debug version:
#    $ make debug
# 4. To clean up:
#    $ make clean
#

# --- Compiler Setup ---
FC = /opt/intel/oneapi/compiler/2025.2/bin/ifx

# --- Target Files ---
# The name of main program
TARGET_NAME = discrete_state
# Directory for source files
SRC_DIR = SourceFiles
# Directory for compiled .o and .mod files
OBJ_DIR = build
# The final executable file
TARGET = $(TARGET_NAME)

# --- Compilation Flags ---

# Standard release flags:
# -O3          	... High optimization level
# -fopenmp     	... Enable OpenMP parallelization
# -mkl=parallel ... Automatically find and link the parallel MKL library
FFLAGS = -O3 -fopenmp -qmkl=parallel -I$(SRC_DIR) -I$(OBJ_DIR)

# Linker flags (MKL and OpenMP are already handled by FFLAGS)
LDFLAGS = -fopenmp -qmkl=parallel

# Debug flags:
# -g           ... Include debugging symbols
# -O0          ... Disable optimization (for easier debugging)
# -traceback   ... Print where the program crashed
# -check all   ... Enable all runtime checks (array bounds, etc.) - very slow!
DEBUG_FLAGS = -g -O0 -traceback -check all -I$(SRC_DIR) -I$(OBJ_DIR)

# Debug linker flags 
DEBUG_LDFLAGS = -fopenmp -qmkl=parallel -g

# --- Automatic File Detection ---
SRCS = $(wildcard $(SRC_DIR)/*.f90)
OBJS = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRCS))

# --- Build Rules ---
all: $(TARGET)

$(TARGET): $(OBJS)
	@echo "==> Linking executable $(TARGET)..."
	$(FC) $(LDFLAGS) -o $(TARGET) $(OBJS)
	@echo "==> Build complete: $(TARGET)"


# --- EXPLICIT DEPENDENCIES ---
build/Math.o: build/Parameters.o
build/Green.o: build/Parameters.o  build/Math.o build/whittaker_w.o build/coulcc.o
build/DiscreteState.o: build/Parameters.o build/Math.o build/Green.o build/whittaker_w.o build/coulcc.o
build/RadialShooting.o: build/Parameters.o build/Math.o build/whittaker_w.o 
build/PHPBound.o: build/Parameters.o build/Math.o build/Green.o build/whittaker_w.o
build/GDTrans.o: build/Parameters.o build/Math.o

build/MAIN.o: build/Parameters.o build/whittaker_w.o build/coulcc.o build/DiscreteState.o build/Green.o build/GDTrans.o build/RadialShooting.o build/PHPBound.o


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(OBJ_DIR)
	@echo "==> Compiling $<..."
	$(FC) $(FFLAGS) -c $< -o $@ -module $(OBJ_DIR)


debug:
	@echo "==> Building DEBUG version..."
	$(MAKE) all FFLAGS="$(DEBUG_FLAGS)" LDFLAGS="$(DEBUG_LDFLAGS)" TARGET_NAME="debug_$(TARGET_NAME)"
	@echo "==> Build complete: debug_$(TARGET_NAME)"


clean:
	@echo "==> Cleaning build files..."
	rm -rf $(OBJ_DIR) $(TARGET) debug_$(TARGET_NAME)
	@echo "==> Cleaning stray .mod files from Source..."
	rm -f $(SRC_DIR)/*.mod