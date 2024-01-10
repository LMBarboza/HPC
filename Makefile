CC = g++

CFLAGS = -mavx -fopenmp

SRC_DIR = src

SRCS = $(SRC_DIR)/main.cpp $(SRC_DIR)/matrix.cpp $(SRC_DIR)/matrixAVX.cpp $(SRC_DIR)/matrixMP.cpp

OUTPUT_DIR = bin

TARGET = $(OUTPUT_DIR)/avx_mp.out

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET)

clean:
	rm -f $(TARGET)
