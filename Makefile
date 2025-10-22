# Compiler settings
CXX = g++
CXXFLAGS = -std=c++26 -Wall -Wextra -O2
INCLUDES = -Ilib

# Directories
SRC_DIR = src
LIB_DIR = lib
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj

# Target executable
TARGET = $(BUILD_DIR)/SPH_Simulation

# Source files
MAIN_SRC = $(SRC_DIR)/main.cpp
LIB_SRCS = $(wildcard $(LIB_DIR)/*.cpp)

# Object files
MAIN_OBJ = $(OBJ_DIR)/main.o
LIB_OBJS = $(patsubst $(LIB_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(LIB_SRCS))

# All object files
OBJS = $(MAIN_OBJ) $(LIB_OBJS)

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJS)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)
	@echo "Build complete: $(TARGET)"

# Compile main.cpp
$(MAIN_OBJ): $(MAIN_SRC)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Compile library files
$(OBJ_DIR)/%.o: $(LIB_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR)
	@echo "Clean complete"

# Run the program
run: $(TARGET)
	./$(TARGET)

# Rebuild everything
rebuild: clean all

# Show variables (for debugging the Makefile)
debug:
	@echo "MAIN_SRC: $(MAIN_SRC)"
	@echo "LIB_SRCS: $(LIB_SRCS)"
	@echo "OBJS: $(OBJS)"

.PHONY: all clean run rebuild debug
