#!/bin/bash

cp ../../commons_portable.* .
g++ mainJ4.C commons_portable.cpp graph_lib.C -O3 -o main
./main
rm main commons_portable.*
