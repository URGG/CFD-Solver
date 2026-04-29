CC = gcc

TARGET_NAME := program

RAYLIB_INC := libraries/raylib/include/
MLIB_INC := libraries/mlib/include/
GSL_INC := libraries/gsl/include/

RAYLIB_LIB := libraries/raylib/liblnx/
GSL_LIB := libraries/gsl/liblnx/

LIBS = -lraylib -lgsl -lgslcblas -lGL -lX11 -lm -lpthread -ldl -lrt

ifeq ($(OS),Windows_NT)
    TARGET_NAME := $(addsuffix .exe,$(TARGET_NAME))
    RAYLIB_LIB := libraries/raylib/libwin/
    GSL_LIB := libraries/gsl/libwin/
    LIBS := -lraylib -lgsl -lgslcblas -lopengl32 -lgdi32 -lwinmm
endif

CFLAGS = -Wall -Wextra -g -I$(RAYLIB_INC) -I$(MLIB_INC) -I$(GSL_INC)
DBGFLAGS := -g
COBJFLAGS := $(CFLAGS) -c
LDFLAGS = -L$(RAYLIB_LIB) -L$(GSL_LIB) $(LIBS)

BIN_PATH := bin
OBJ_PATH := obj
SRC_PATH := src
DBG_PATH := debug

TARGET := $(BIN_PATH)/$(TARGET_NAME)
TARGET_DEBUG := $(DBG_PATH)/$(TARGET_NAME)

SRC := $(foreach x, $(SRC_PATH), $(wildcard $(addprefix $(x)/*,.c*)))
OBJ := $(addprefix $(OBJ_PATH)/, $(addsuffix .o, $(notdir $(basename $(SRC)))))
OBJ_DEBUG := $(addprefix $(DBG_PATH)/, $(addsuffix .o, $(notdir $(basename $(SRC)))))

DISTCLEAN_LIST := $(OBJ) \
                  $(OBJ_DEBUG)
CLEAN_LIST := $(TARGET) \
            $(TARGET_DEBUG) \
            $(DISTCLEAN_LIST)

default: makedir all

$(TARGET): $(OBJ)
	$(CC) -o $@ $(OBJ) $(CFLAGS) $(LDFLAGS)

# Removed LDFLAGS from the object compilation steps
$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c*
	$(CC) $(COBJFLAGS) -o $@ $<

$(DBG_PATH)/%.o: $(SRC_PATH)/%.c*
	$(CC) $(COBJFLAGS) $(DBGFLAGS) -o $@ $<

# Moved LDFLAGS to the end of the line for correct linking order
$(TARGET_DEBUG): $(OBJ_DEBUG)
	$(CC) $(CFLAGS) $(DBGFLAGS) $(OBJ_DEBUG) -o $@ $(LDFLAGS)

.PHONY: makedir
makedir:
	@mkdir -p $(BIN_PATH) $(OBJ_PATH) $(DBG_PATH)

.PHONY: all
all: $(TARGET)

.PHONY: debug
debug: $(TARGET_DEBUG)

.PHONY: clean
clean:
	@echo CLEAN $(CLEAN_LIST)
	@rm -f $(CLEAN_LIST)

.PHONY: distclean
distclean:
	@echo CLEAN $(DISTCLEAN_LIST)
	@rm -f $(DISTCLEAN_LIST)

run:
	make clean
	make
	./$(BIN_PATH)/$(TARGET_NAME)