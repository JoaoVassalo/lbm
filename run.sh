#!/bin/bash

set -e

NVCC=/usr/local/cuda/bin/nvcc

SRC_DIR=src
BUILD_DIR=build
BIN=debug

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