@echo off

cls
gcc -shared -fpic -O3 -ffast-math -o tridiag.dll cochlea_utils.c 
