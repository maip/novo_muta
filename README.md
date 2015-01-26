novo_muta
=========

@author Melissa Ip

An implementation of probabilistic methods for detecting de novo mutations from nuclear family genome data.

See Cartwright et al.: Family-Based Method for Capturing De Novo Mutations:
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3728889/

There are currently two versions of the trio model. One uses customized Dirichlet-multinomial approximations and has a simulation. The trio model based on the infinite sites model uses simpler multinomial approximations and does not have a simulation.

#Required libraries

This project uses C++11, Eigen, GNU Scientific Library, BamTools, and Boost. It compiles using Clang and CMake, configured specifically for use with herschel (Cartwright lab). You may need to modify ```CMakeLists.txt``` to make this project compile for your system.

Download [Eigen](http://eigen.tuxfamily.org/) and put the source code in the project directory.

Download [GSL](http://www.gnu.org/software/gsl/).

Download [BamTools](https://github.com/pezmaster31/bamtools). Keep a local copy of the ```utils``` folder in the project directory as the file ```utils/bamtools_pileup_engine.h``` is required.

Download [Boost](http://www.boost.org/users/download/). Be sure to include and build the unit test framework.

#Compilation

```mkdir build``` makes a build directory.

```cd build``` navigates to the new build directory.

```cmake ..``` reads ```CMakeLists.txt``` and makes the build files. ```FindGSL.cmake``` and ```FindBamTools.cmake``` are provided in the ```Modules``` folder.

```make``` builds the CMake files.

To execute a file, use ```./<filename>```.

#Testing

Navigate to the build directory you created earlier to compile the project and type ```make test``` to run all Boost tests in the test folder. If there are any changes to ```CMakeLists.txt```, the command ```cmake ..``` must precede ```make test``` in order to update.

Individual tests can be executed using ```./<filename>```.
