# CLOVERS Makefile

CXX := g++
CXXFLAGS ?= -DZLIB -fopenmp -mavx -mfma -static -O3
LDFLAGS ?= -lz

SRC_DIR := src
INC_DIR := include
BUILD_DIR := build

# Object files (auto-derived from src/*.cpp, preserving filename casing)
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(SRC_FILES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Target executables (per-platform filenames)
ifeq ($(OS),Windows_NT)
    TARGET_FILES := clovers.exe labeler.exe
    EXT := .exe
    CLEAN_CMD := powershell -Command " \
        Remove-Item -Force -ErrorAction SilentlyContinue clovers.exe, labeler.exe; \
        Remove-Item -Recurse -Force -ErrorAction SilentlyContinue '$(BUILD_DIR)'"
else
    TARGET_FILES := clovers labeler
    EXT :=
    CLEAN_CMD := rm -f $(TARGET_FILES) && rm -rf $(BUILD_DIR)
endif

.PHONY: all clean rebuild

all: $(BUILD_DIR) $(TARGET_FILES)

$(BUILD_DIR):
	mkdir $(BUILD_DIR)

# Compile any .cpp to .o (casing matches source filenames)
$(OBJ_FILES): $(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

# Per-target prerequisites
clovers$(EXT): $(BUILD_DIR)/svm.o $(BUILD_DIR)/Clovers.o
labeler$(EXT): $(BUILD_DIR)/svm.o $(BUILD_DIR)/Labeler.o

# Shared link recipe
$(TARGET_FILES):
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	$(CLEAN_CMD)

rebuild: clean all
