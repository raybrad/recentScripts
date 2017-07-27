#!/bin/sh
awk '/k-points in units of 2pi\/SCALE and weight:/,/k-points in reciprocal lattice and weights:/' OUTCAR > kpoints.dat
