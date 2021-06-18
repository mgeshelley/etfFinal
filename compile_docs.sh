#!/bin/bash

cp docs/Doxyfile . # Get Doxyfile from docs directory
doxygen Doxyfile # Compile Doxygen documentation

cd latex
make # Compile LaTeX Doxygen documentation
cd ..

# Cleanup
rm -r docs/html
rm -r docs/latex
mv html docs
mv latex docs
rm Doxyfile 