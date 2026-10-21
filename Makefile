# CLOVERS Makefile
#
# Build:    make
#           make OPENMP=0            no OpenMP (no libgomp/libomp needed)
#           make NATIVE=1            AVX/FMA/AVX2 (pre-Haswell CPUs crash)
#           make STATIC=1            static binary
# Install:  make install PREFIX=/usr
#
# Overridable: CXX CPPFLAGS CXXFLAGS LDFLAGS LDLIBS PREFIX DESTDIR INSTALL

SRC_DIR     := src
INC_DIR     := include
BUILD_DIR   := build
SCRIPTS_DIR := scripts

CXX      ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -O3
LDFLAGS  ?=
LDLIBS   ?= -lz
INSTALL  ?= install
PREFIX   ?= /usr/local
DESTDIR  ?=

CXXFLAGS += -std=c++17
CPPFLAGS += -DZLIB

OPENMP ?= 1
STATIC ?= 0
NATIVE ?= 0
OPENMP_FLAGS ?= -fopenmp
OPENMP_LIBS  ?= -fopenmp

ifeq ($(OPENMP),1)
    CXXFLAGS += $(OPENMP_FLAGS)
    LDFLAGS  += $(OPENMP_LIBS)
endif

ifeq ($(NATIVE),1)
    CXXFLAGS += -mavx -mfma -mavx2
endif

ifeq ($(OS),Windows_NT)
    EXT       := .exe
    MKDIR_P   := mkdir
    STATIC    ?= 1          # MinGW needs the C++/gomp runtimes bundled
    CLEAN_CMD = powershell -Command " \
        Remove-Item -Force -ErrorAction SilentlyContinue clovers.exe; \
        Remove-Item -Force -ErrorAction SilentlyContinue '$(SCRIPTS_DIR)/labeler.exe'; \
        Remove-Item -Recurse -Force -ErrorAction SilentlyContinue '$(BUILD_DIR)'"
else
    EXT       :=
    MKDIR_P   := mkdir -p
    CLEAN_CMD = rm -f clovers $(SCRIPTS_DIR)/labeler && rm -rf $(BUILD_DIR)
endif

ifeq ($(STATIC),1)
    LDFLAGS += -static
endif

SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(SRC_FILES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
HEADERS   := $(wildcard $(INC_DIR)/*.hpp)

CLOVERS := clovers$(EXT)
LABELER := $(SCRIPTS_DIR)/labeler$(EXT)
TARGETS := $(CLOVERS) $(LABELER)

.PHONY: all clean rebuild install native

all: $(TARGETS)

$(BUILD_DIR) $(SCRIPTS_DIR):
	$(MKDIR_P) $@

$(OBJ_FILES): $(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

$(CLOVERS): $(BUILD_DIR)/svm.o $(BUILD_DIR)/Clovers.o
$(LABELER): $(BUILD_DIR)/svm.o $(BUILD_DIR)/Labeler.o

$(TARGETS): | $(SCRIPTS_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

install: $(CLOVERS)
	$(INSTALL) -d $(DESTDIR)$(PREFIX)/bin
	$(INSTALL) -m 755 $(CLOVERS) $(DESTDIR)$(PREFIX)/bin/

native:
	$(MAKE) NATIVE=1 STATIC=1

clean:
	$(CLEAN_CMD)

rebuild: clean all
