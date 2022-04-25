#!/bin/bash
  echo -e "Size\t\tNo-opt\t\tCode\t\tNo-mem\t\t2x1\t\t2x2\t\t4x1\t\t8x1"
  for n in {10..50..10}; do
     gcc -O0 -DN=$n -oprog2 prog2.c -lm
     ./prog2
  done
