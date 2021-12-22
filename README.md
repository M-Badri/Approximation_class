# Approximation Class


## This library solves the approximation (interpolation) problems in 2D and 3D space. 
* This library solves the approximation (interpolation) problems in 2D and 3D space.
* Although theoritcally has no limitation, it has been desinged for AMR and OverSet Mesh methods maily.

 
## How to use

### Using the source code direvctly in your projects.
* We have added all needed source codes to the library, so hopefly it can be independently used.


### Using Cmake (assuming that cmake library is already installed)
* We have provided a cmake source file for the project that can be employed as ususal. At the command prompt type:
  * `mkdir build`
  * `cd build` 
  * `cmake .. -D CMAKE_Fortran_COMPILER=ifort`
  * **NOTE** The makefile is made in the _build_. You can make the project just by typing `make` in the _build_ directory.
  * Then you can find the executable files in the _build/bin_ directory.
  
<!---  
* Additionally, you can use cmake to build a CodeBlocks prject and then work with CodeBlocks IDE. At the command prompt type:
 	* `cmake . -G "CodeBlocks - Unix Makefiles" -D CMAKE_Fortran_COMPILER=Intel`
 	* After that you can find a _cbp_ file of the CodeBlocks in the main directory and it can be used as usual CB project. 
 	* **NOTE** While working with CB if you get the error `you must select host application to run ...` the fllowing steps may help to get around it:
 		+ By right clicking on the project 
 		+ project > properties > Build targets
 		+ From the _Type_ drop down list select _Console Application_
 		+ In the output filename section write: bin/the_executable_file_name (in my case: bin/example_1_2d)
 -->

	


