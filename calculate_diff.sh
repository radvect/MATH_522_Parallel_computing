#!/bin/bash


MATRIX1="grid_parallel.txt" 
MATRIX2="grid_nonparallel.txt"

MAX_NORM=$(paste "$MATRIX1" "$MATRIX2" | awk '{
    max = 0;
    for (i=2; i<=NF/2-1; i++) {
        diff = ($i - $(i + NF/2));
        if (diff < 0) diff = -diff;
        if (diff > max) max = diff;
    }
} END { print max }')

echo "L_\inf norm: $MAX_NORM"