#!/bin/bash

mkdir 2D2Gdiff
cp -r include 2D2Gdiff
cp -r src 2D2Gdiff
cp -r 2D_1G.inp 2D2Gdiff
cp -r Makefile 2D2Gdiff
cp -r make_dependencies.pl 2D2Gdiff
cp -r plotflux.py 2D2Gdiff
cp -r matplot.py 2D2Gdiff

tar -zcvf 2D2Gdiff.tar.gz 2D2Gdiff
rm -rf 2D2Gdiff
