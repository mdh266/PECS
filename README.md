## Introduction

This software is designed to simulate the dynamics of the reactive interface between 
a semiconductor and electrolyte of a photoelectrochemical solar cell. 

For **much more** background on the model and algorithms used in this project please see the 
 <a href="http://mdh266.github.io/PECS/">documentation page</a>.

## Requirements

The requirements for this software are <a href="dealii.org">deal.ii</a> library version 
8.3.0 or higher, and <a href="https://cmake.org/">CMake</a> version 2.8 or higher. The code 
will automatically run in parallel using the 
<a href="https://www.threadingbuildingblocks.org/">Thread Building Blocks</a>. 

## Using

You need to obtain and install a copy of the deal.ii library version 8.3.0 or higher. 
After downloading and installing the deal.II library. cd into the PECS directory.

To generate a make file run to compile the source code type:

	cmake . -DDEAL_II_DIR="path to deal.II library"	

On a mac, if you downloaded the binaries of deal.ii library instead run:

	cmake .

Once this complete you can type:

	make release	

to compile the code.

To run the code type:

	./main	

Material and design choices can be chosen by the user through the input file,

	input_file.prm	

See the documentation for information on the input parameter values.

The resulting outputs are in <a href="http://www.vtk.org/">VTK</a> format and can viewed using 
Paraview<a href="http://www.paraview.org/">Paraview</a>.



