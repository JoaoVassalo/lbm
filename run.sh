#!/bin/bash

set -e

NVCC=/usr/local/cuda-13.1/bin/nvcc

SRC_DIR=src
BUILD_DIR=build
BIN=debug

echo "Limpando..."

# Remove o executável antigo
rm -f ${BUILD_DIR}/${BIN}

# Remove os arquivos VTK antigos
rm -rf results/plot/vtk/*

mkdir -p ${BUILD_DIR}

echo "Compilando..."

$NVCC \
    -std=c++20 \
    -g -G \
    -rdc=true \
    -I${SRC_DIR} \
    $(find ${SRC_DIR} -name "*.cu") \
    -o ${BUILD_DIR}/${BIN}

echo "Executando..."

./${BUILD_DIR}/${BIN}