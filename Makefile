TARGET ?= Phase-field
FOLDER ?= Results

BUILD_DIR := ./build
INCLUDE_DIR := -I./include -I~/.local/include
SRC_DIR := ./src
RESULTS_DIR := ./results
PYTHON_DIR := ./python
PARAM_PATH := info/parameters.txt

CONFIGS_DIR := ./comp_configurations
CONFIG_TARGET := ./config/config.json

COMPILER = nvcc
COMPILER_FLAGS = -std=c++17 $(INCLUDE_DIR) -O3 --expt-relaxed-constexpr --extended-lambda -rdc=true $(TNL_FLAGS)

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
	echo "" >> $(RESULTS_DIR)/$(FOLDER)/$(PARAM_PATH)
	echo -n "Run date: " >> $(RESULTS_DIR)/$(FOLDER)/$(PARAM_PATH)
	date "+%d.%m.%Y %H:%M:%S" >> $(RESULTS_DIR)/$(FOLDER)/$(PARAM_PATH)
	echo -n "Git commit: " >> $(RESULTS_DIR)/$(FOLDER)/$(PARAM_PATH)
	git rev-parse --short HEAD >> $(RESULTS_DIR)/$(FOLDER)/$(PARAM_PATH)
	python3 $(PYTHON_DIR)/basic_plot.py --name $(FOLDER) --type both
	dir="$(RESULTS_DIR)/$(FOLDER)/calculations"; \
	if [ -d "$$dir" ]; then \
		files=($$(ls -1 "$$dir")); \
		count=$${#files[@]}; \
		if [ $$count -gt 2 ]; then \
			to_delete=$$(printf "%s\n" "$${files[@]:1:$$(($$count-2))}"); \
			for f in $$to_delete; do rm -rf "$$dir/$$f"; done; \
		fi; \
	fi; \

.PHONY: run_all
run_all:
	@for cfg in $(CONFIGS_DIR)/*.json; do \
		name=$$(basename $$cfg .json); \
		echo "Spouštím výpočet pro konfiguraci: $$name"; \
		echo "========================================"; \
		cp $$cfg $(CONFIG_TARGET); \
		$(MAKE) run FOLDER=$$name; \
	done

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)/*
