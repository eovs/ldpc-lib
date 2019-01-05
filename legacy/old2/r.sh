#!/bin/bash

gcc -O3 -o main mainJ4.C graph_mem.C graph_lib.C
./main
rm main
