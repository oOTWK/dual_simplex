CC = gcc

CFLAGS = -Wall -O3 -std=gnu99

BUILD_DIR = build

OBJ = $(BUILD_DIR)/main.o $(BUILD_DIR)/dual_simplex.o

$(BUILD_DIR)/bin/dual_simplex: $(OBJ)  	
	@ echo Linking Binary: $@
	@ mkdir -p $(BUILD_DIR)/bin
	@ $(CC) $(CFLAGS) $^ -lm -o $@

$(BUILD_DIR)/main.o: main.c
	@ echo Compiling: $@
	@ mkdir -p $(BUILD_DIR)
	@ $(CC) $(CFLAGS) $< -c -o $@

$(BUILD_DIR)/dual_simplex.o: dual_simplex.c
	@ echo Compiling: $@
	@ mkdir -p $(BUILD_DIR)
	@ $(CC) $(CFLAGS) $< -lm -c -o $@


.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)