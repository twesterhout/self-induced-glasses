#!/bin/bash

for k in $(seq 0 0); do
  for i in $(seq 1 60); do
    /vol/tcm40/westerhout_tom/self-induced-glasses/dist-newstyle/build/x86_64-linux/ghc-8.10.7/self-induced-glasses-0.0.0.0/x/self-induced-glasses/opt/build/self-induced-glasses/self-induced-glasses \
      +RTS -N1 -RTS -n 25 -λ 1.5 --seed $((47607 + i + 60 * k)) --sweep-size 40000 \
      --step 0.2,100,1000 \
      --step 0.6,100,1000 \
      --step 1.0,100,1000 \
      --step 1.4,100,1000 \
      --step 1.8,100,1000 \
      --step 2.2,100,1000 \
      --step 2.6,100,1000 \
      --step 3.0,100,1000 &
  done
  wait
done
