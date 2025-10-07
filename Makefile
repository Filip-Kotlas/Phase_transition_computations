TARGET ?= Phase-field
FOLDER ?= Results

BUILD_DIR := ./build
INCLUDE_DIR := -I./include -I~/.local/include
SRC_DIR := ./src
RESULTS_DIR := ./results
PYTHON_DIR := ./python

COMPILER = nvcc
COMPILER_FLAGS = -std=c++17 $(INCLUDE_DIR) -O3 --expt-relaxed-constexpr --extended-lambda $(TNL_FLAGS)

# Potlačení zbytečných varování
# -diag-suppress 20012 = "annotation ignored" (ArrayView apod.)
# -diag-suppress 177 = "variable declared but never referenced"
# -diag-suppress 186 = "pointless comparison"
# -diag-suppress 297 = "declared but never referenced"
# -Xcompiler -w = potlačí běžná CPU varování (gcc/clang)
SUPPRESS_WARNINGS := -Xcudafe --diag_suppress=20012 -Xcudafe --diag_suppress=177 \
                     -Xcudafe --diag_suppress=186 -Xcudafe --diag_suppress=297 \
                     -Xcompiler -w

SRC := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*.cu) $(wildcard $(SRC_DIR)/*.c)
OBJ := $(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(SRC:.cpp=.o))
OBJ := $(OBJ:.cu=.o)
OBJ := $(OBJ:.c=.o)

# pravidla pro překlad
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(COMPILER) $(COMPILER_FLAGS) $(SUPPRESS_WARNINGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu
	@mkdir -p $(BUILD_DIR)
	$(COMPILER) $(COMPILER_FLAGS) $(SUPPRESS_WARNINGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(COMPILER) $(COMPILER_FLAGS) $(SUPPRESS_WARNINGS) -x c -c $< -o $@

# linkování
$(BUILD_DIR)/$(TARGET): $(OBJ)
	$(COMPILER) $(COMPILER_FLAGS) $(SUPPRESS_WARNINGS) $(OBJ) -o $@ $(TNL_LDFLAGS)

.PHONY: run
run: $(BUILD_DIR)/$(TARGET)
	$(BUILD_DIR)/$(TARGET) $(RESULTS_DIR)/$(FOLDER)
#python3 $(PYTHON_DIR)/basic_plot.py --name $(FOLDER) --type both

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)/*
