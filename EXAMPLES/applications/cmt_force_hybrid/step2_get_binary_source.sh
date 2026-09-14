#!/bin/bash
gfortran write_binary_source.f90 -o write_source

./write_source
\rm write_source