**Introduction**

This software is designed to simulate the dynamics of the reactive interface between a semiconductor and electrolyte. The interface of between the semiconductor and electrolyte make up a "half cell" of a photoelectrochemical cell which can use solar energy to convert water into hydrogen fuel.

**Background and Documentation**
For more information on this project please see the documenation page: <a href="http://mdh266.github.io/PECS/">http://michael-harmon.com/PECS/</a>

**Requirements**

The requirements for this software are <a href="dealii.org">deal.ii</a> library version 8.3.0 or higher, and <a href="https://cmake.org/">CMake</a> version 2.8 or higher. The code will automatically run in parallel using the Thread Building Blocks. See deal.ii's explanation on parallel computing with shared memory to see why this necessary and how it works. To make this documentation on your local machine you need doxygen.

**Usage**

You need to obtain and install a copy of the dealii deal.ii library version 8.3.0 or higher. After downloading and installing the deal.II library. cd into the PECS directory.

To generate a make file run to compile the source code type:

*cmake . -DDEAL_II_DIR="path to deal.II library"*

On a mac, if you downloaded the binaries of deal.ii library instead run:

*cmake .*

Once this complete you can type:

*make release*

to compile the code.

To run the code type:

*./main*

Material and design choices can be chosen by the user through the input file input_file.prm.

The resulting outputs are in VTK format and can viewed using Paraview.


