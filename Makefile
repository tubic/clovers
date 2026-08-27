# CLOVERS Makefile

CXX := g++
CXXFLAGS ?= -DZLIB -fopenmp -mavx -mfma -static -O3
LDFLAGS ?= -lz

SRC_DIR := src
INC_DIR := include
BUILD_DIR := build
SCRIPTS_DIR := scripts

# Object files (auto-derived from src/*.cpp, preserving filename casing)
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(SRC_FILES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Per-platform filenames
ifeq ($(OS),Windows_NT)
    EXT := .exe
    CLEAN_CMD := powershell -Command " \
        Remove-Item -Force -ErrorAction SilentlyContinue clovers.exe; \
        Remove-Item -Force -ErrorAction SilentlyContinue '$(SCRIPTS_DIR)/labeler.exe'; \
        Remove-Item -Recurse -Force -ErrorAction SilentlyContinue '$(BUILD_DIR)'"
else
    EXT :=
    CLEAN_CMD := rm -f clovers $(SCRIPTS_DIR)/labeler && rm -rf $(BUILD_DIR)
endif

# Target executables (paths to the final binaries)
CLOVERS := clovers$(EXT)
LABELER := $(SCRIPTS_DIR)/labeler$(EXT)
TARGET_FILES := $(CLOVERS) $(LABELER)

.PHONY: all clean rebuild

all: $(BUILD_DIR) $(TARGET_FILES)

$(BUILD_DIR):
	mkdir $(BUILD_DIR)

# Compile any .cpp to .o (casing matches source filenames)
$(OBJ_FILES): $(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

# Per-target prerequisites
$(CLOVERS): $(BUILD_DIR)/svm.o $(BUILD_DIR)/Clovers.o
$(LABELER): $(BUILD_DIR)/svm.o $(BUILD_DIR)/Labeler.o

# Shared link recipe (output path is the target name)
$(TARGET_FILES):
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	$(CLEAN_CMD)

rebuild: clean all
