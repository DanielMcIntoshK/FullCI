SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin
INC_DIR := include
EXE := $(BIN_DIR)/DanielsFullCI
SRC := $(wildcard $(SRC_DIR)/*.cxx)
OBJ := $(SRC:$(SRC_DIR)/%.cxx=$(OBJ_DIR)/%.o)
CPPFLAGS := -I$(INC_DIR)

all: $(EXE)

CC := g++

.PHONY: all

$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) -O2 $^ -o $@
$(BIN_DIR):
	$(CC) -p $@
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx | $(OBJ_DIR) 
	$(CC) $(CPPFLAGS) -c $< -o $@
$(OBJ_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR)
