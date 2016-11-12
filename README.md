This is a code for simulating the dynamics of a photoelectrochemical solar cell in two dimensions. 



**Dependencies**

It requires the 
[CMake (v 2.8.8) ](http://cmake.org), the [deal.II library (v 8.4)](http://dealii.org) 
to run and [Paraview (v 5.0)](http://www.paraview.org) for visualization.  

**Installation**

After downloading and installing the deal.II library.  cd into pecs directory.  
To generate a make file run: 

*cmake . -DDEAL_II_DIR="path to deal.II library"*

Once this complete you can type:

*make release*

to compile the code. Note: if just use:

*make* 

then there will be issues with the printing, if you want to use that version,
 you need to uncomment the commented part of the printing call.


**Running The Code**

To run the code type:

*./main*

To clean up the runtime output files type:

*make runclean*


**Documentation**

For more information cd into "/Documentation" and run "doxygen dox" 
to generate thorough html documentation.


**Questions**

Email Mike: mdh266@gmail.com
