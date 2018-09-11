#!/bin/bash

clear
gcc -shared -fpic -O3 -ffast-math -o tridiag.so cochlea_utils.c 
