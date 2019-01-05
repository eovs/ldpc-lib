#!/bin/bash

g++ -O3 commons_portable.cpp graph_lib.C graph_mem.C mainJ4.C -o mainJ4
time ./mainJ4 input32_16.cfg  out.dat
rm mainJ4
