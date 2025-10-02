TARGET ?= Phase-field
FOLDER ?= Results

BUILD_DIR := ./build
INCLUDE_DIR := ./include
SRC_DIR := ./src
RESULTS_DIR := ./results
PYTHON_DIR := ./python

CXX = g++
CXXFLAGS = -std=c++17 -I$(INCLUDE_DIR)
WARNINGS := -Wall -Wextra

SRC := $(wildcard $(SRC_DIR)/*.cpp)
INCLUDE := $(wildcard $(INCLUDE_DIR)/*.hpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# pravidlo pro překlad .cpp -> .o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# linkování
$(BUILD_DIR)/$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: run
run: $(BUILD_DIR)/$(TARGET)
	$(BUILD_DIR)/$(TARGET) $(RESULTS_DIR)/$(FOLDER)
	python3 $(PYTHON_DIR)/basic_plot.py --name $(FOLDER) --type both

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)/*
