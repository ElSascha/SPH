# Compiler settings
CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra -O2
INCLUDES = -Ilib -Ilib/HighFive/include
LDFLAGS = -lhdf5_cpp -lhdf5

# Directories
SRC_DIR = src
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin

# Target executable
TARGET = $(BIN_DIR)/SPH_Simulation

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)
	@echo "✅ Build complete: $(TARGET)"

# Compile source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR)
	@echo "🧹 Clean complete"

# Run the program
run: $(TARGET)
	./$(TARGET)

# Rebuild everything
rebuild: clean all

debug:
	@echo "SRCS: $(SRCS)"
	@echo "OBJS: $(OBJS)"