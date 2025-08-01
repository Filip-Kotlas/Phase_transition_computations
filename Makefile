TARGET ?= Phase-field.exe
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

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@if not exist "$(BUILD_DIR)" mkdir $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY:run
run: $(BUILD_DIR)/$(TARGET)
	$(BUILD_DIR)/$(TARGET) $(RESULTS_DIR)/$(FOLDER)
	python $(PYTHON_DIR)/basic_plot.py --name $(FOLDER) --type both

.PHONY:clean
clean:
	if exist "$(BUILD_DIR)\*" del /Q /S "$(BUILD_DIR)\*"