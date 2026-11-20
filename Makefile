TARGET ?= Phase-field
FOLDER ?= Results

BUILD_ROOT := ./build
# Compute capability of the GPU on the build machine (e.g. 86). Objects are kept in a
# separate folder per GPU type, so a shared home directory cannot mix up architectures.
GPU_TAG := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -n 1 | tr -d '. \r')
BUILD_DIR := $(BUILD_ROOT)/sm_$(or $(GPU_TAG),generic)
INCLUDE_DIR := -I./include -I~/.local/include
SRC_DIR := ./src
RESULTS_DIR := ./results
PYTHON_DIR := ./python
PARAM_PATH := info/setup/parameters.txt

DELETE_CALCS ?= 0

CONFIGS_DIR := ./comp_configurations
CONFIG_TARGET := ./config/config.json
# Commit recorded in parameters.txt. On the cluster there is no .git folder, so cluster.sh passes it in.
GIT_COMMIT ?= $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
# Config file the program reads (make run CONFIG=path/to/other.json)
CONFIG ?= $(CONFIG_TARGET)

# Compile for the GPU that is present on the build machine (needs a reasonably recent nvcc).
# To force one architecture: make CUDA_ARCH="-gencode arch=compute_86,code=sm_86"
CUDA_ARCH ?= -arch=native

COMPILER = nvcc
COMPILER_FLAGS = -std=c++17 $(INCLUDE_DIR) -O3 --expt-relaxed-constexpr --extended-lambda -rdc=true $(CUDA_ARCH)

# Suppress unnecessary warnings
# -diag-suppress 20012 = "annotation ignored" (ArrayView etc.)
# -diag-suppress 177 = "variable declared but never referenced"
# -diag-suppress 186 = "pointless comparison"
# -diag-suppress 297 = "declared but never referenced"
# -Xcompiler -w = suppresses common CPU warnings (gcc/clang)
SUPPRESS_WARNINGS := -Xcudafe --diag_suppress=20012 -Xcudafe --diag_suppress=177 \
                     -Xcudafe --diag_suppress=186 -Xcudafe --diag_suppress=297 \
                     -Xcompiler -w

# All sources, also in the subfolders of src/ (config/, models/, ...). The objects keep the same
# subfolders in the build folder, e.g. build/sm_75/models/Model.o.
SRC := $(sort $(shell find $(SRC_DIR) -name '*.cpp' -o -name '*.cu' -o -name '*.c'))
OBJ := $(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(SRC:.cpp=.o))
OBJ := $(OBJ:.cu=.o)
OBJ := $(OBJ:.c=.o)

# Header dependencies: while compiling, nvcc also writes build/sm_XX/<name>.d listing the
# headers the object includes (-MMD: only the project's headers, not system ones such as TNL
# in /usr/include). These files are included at the end of this Makefile, so changing a
# header recompiles exactly the objects that include it. -MP adds an empty rule for each
# header, so deleting or renaming a header does not break the build.
DEPFLAGS = -MMD -MP -MF $(@:.o=.d)

# compilation rules
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(COMPILER) $(COMPILER_FLAGS) $(SUPPRESS_WARNINGS) $(DEPFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu
	@mkdir -p $(@D)
	$(COMPILER) $(COMPILER_FLAGS) $(SUPPRESS_WARNINGS) $(DEPFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(COMPILER) $(COMPILER_FLAGS) $(SUPPRESS_WARNINGS) $(DEPFLAGS) -x c -c $< -o $@

# linking
$(BUILD_DIR)/$(TARGET): $(OBJ)
	$(COMPILER) $(COMPILER_FLAGS) $(SUPPRESS_WARNINGS) $(OBJ) -o $@ $(TNL_LDFLAGS)

# The binary `make run` runs. On the cluster, run.sbatch passes the one from the build cache
# (slurm/cached_build.sh); make has no rule for a file outside build/, so it does not rebuild it.
BINARY ?= $(BUILD_DIR)/$(TARGET)

.PHONY: run
run: $(BINARY)
	@# Exit code 2 (EXIT_RUN_FAILED in main.cu): the solution broke down (NaN or inf, or the solver got stuck).
	@# The frames up to the failure are still plotted (--force: basic_plot.py refuses an incomplete
	@# run without it), then make fails with that code, so the frames are not deleted (DELETE_CALCS)
	@# and the job counts as failed.
	status=0; $(BINARY) $(RESULTS_DIR)/$(FOLDER) $(CONFIG) || status=$$?; \
	if [ $$status -ne 0 ] && [ $$status -ne 2 ]; then exit $$status; fi; \
	force=""; if [ $$status -eq 2 ]; then force=--force; fi; \
	echo "Git commit: $(GIT_COMMIT)" >> $(RESULTS_DIR)/$(FOLDER)/$(PARAM_PATH); \
	python3 $(PYTHON_DIR)/basic_plot.py --name $(FOLDER) --type both $$force || echo "WARNING: plotting failed, continuing"; \
	exit $$status
	@if [ "$(DELETE_CALCS)" = "1" ]; then \
		dir="$(RESULTS_DIR)/$(FOLDER)/calculations"; \
		if [ -d "$$dir" ]; then \
			files=($$(ls -1 "$$dir")); \
			count=$${#files[@]}; \
			if [ $$count -gt 2 ]; then \
				to_delete=$$(printf "%s\n" "$${files[@]:1:$$(($$count-2))}"); \
				for f in $$to_delete; do rm -rf "$$dir/$$f"; done; \
			fi; \
		fi; \
	fi

.PHONY: run_all
run_all:
	@for cfg in $(CONFIGS_DIR)/*.json; do \
		name=$$(basename $$cfg .json); \
		echo "Running computation for configuration: $$name"; \
		echo "========================================"; \
		cp $$cfg $(CONFIG_TARGET); \
		$(MAKE) run FOLDER=$$name; \
	done

.PHONY: clean
clean:
	rm -rf $(BUILD_ROOT)/*

# Generated header dependencies (see DEPFLAGS). Kept at the end so that they cannot become the
# default goal; the leading "-" ignores the files that do not exist yet (first build).
-include $(OBJ:.o=.d)
