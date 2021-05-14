#!/bin/bash

nvcc -o FWAlgoCUDA FloydWarshall.cu -ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29910\bin\Hostx64\x64"

TIMEFORMAT="%R"
for iter in 1 2 3 4 5 6 7 8 9 10 11
do
  for N in 10 20 40 80 160 320 640 1280 2560 5120
  do
    sleep 0.1
    time ./FWAlgoCUDA $N
  done
done
