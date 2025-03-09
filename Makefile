TARGET_NAME:= Phase-field.exe
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

$(BUILD_DIR)./$(TARGET_NAME): $(SRC) $(INCLUDE)
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)./$(TARGET_NAME) $(SRC)

.PHONY:run
run: $(BUILD_DIR)./$(TARGET_NAME)
	$(BUILD_DIR)./$(TARGET_NAME) $(RESULTS_DIR)/$(FOLDER)
	python $(PYTHON_DIR)/basic_plot.py $(FOLDER)

.PHONY:clean
clean:
	if exist "$(BUILD_DIR)\*" del /Q /S "$(BUILD_DIR)\*"