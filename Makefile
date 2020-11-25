CC = gcc
NVCC = nvcc
FLAGS = -fopenmp -O3 -ffast-math -march=native -fPIC
NVCCFLAGS = -O3 -arch sm_75

all: clean transpose transpose_gpu

clean:
	touch transpose; rm transpose; touch transpose_gpu; rm transpose_gpu

transpose:
	$(CC) $(FLAGS) transpose.c -o transpose

transpose_gpu:
	$(NVCC) $(NVCCFLAGS) transpose_gpu.cu -o transpose_gpu
