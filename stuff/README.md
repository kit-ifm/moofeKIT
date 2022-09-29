# pyfemkit (python finite element kit founded at the karlsruhe intitute of technology)

**main targets of the code**
- python based finite element code 
    - modular
    - minimalistic
    - lightweight
    - object oriented
    - parallel
    - symbolic computations if possible
- keep it simple and stupid (KISS) makes possible the use in advanced courses


**code structure**
- preprocessing
    - preprocessor
    - mesh generators
    - ABAQUS-mesh (.inp-files) interface
    - startscripts
- solver
    - time loop
    - Newton's method
    - fe-assembly
    - direct solver
- postprocessing
    - element plotter
    - diagram plotter (energy, momentum, ...)
    - stress/temperature plotter or Matlab interface


**todo folder (drop tasks with assignments if possible)**
- MATLAB 
    - [x] fe code of finite element course (including rudimentary parallelization) _@Marlon_
    - [x] symbolic toolbox example _@Marlon_
    - [x] add some advanced (mixed) fe code in MATLAB _@Marlon
- python 
    - [x] fe code of finite element course (including rudimentary parallelization) _@Felix_
    - [x] symbolic sympy example _@Felix_
    - [x] test and optimize parallelization _@Felix/Marlon_
    - [x] translate advanced fe (mixed fe) code (start with Aufgabe5b.m which needs elementSRReduced.m and nachlaufSRReduced.m afterwards translate Aufgabe5a.m which needs elementSRFull.m and nachlaufSRFull.m) in MATLAB _@Tim_
- [x] new name for the repo: pyfemkit _@Marlon_
- [x] structure for pyfemkit _@Marlon_
- [x] object oriented architecture _@Marlon_
- [ ] translate ABAQUS interface from MATLAB to python _@Marlon/Felix/Tim_


**literature and links**
- package manager https://www.anaconda.com/products/individual
- integrated development environment (IDE) https://www.spyder-ide.org/
- https://www.w3schools.com/python/numpy_intro.asp
- https://realpython.com/matlab-vs-python/


**Installing via conda-environment e.g. for arch linux**
- yaourt -S anaconda
- optional 
    - source /opt/anaconda/bin/activate
    - conda init
- conda create -n spyder-env spyder=4 numpy scipy pandas matplotlib sympy cython
- sparse direct solver library scikit-umfpack (optional) for better performance
    1. conda install -c conda-forge scikit-umfpack
    2. conda install -c conda-forge/label/gcc7 scikit-umfpack
    3. conda install -c conda-forge/label/cf201901 scikit-umfpack
    4. conda install -c conda-forge/label/cf202003 scikit-umfpack
- conda activate spyder-env (if not already)
- addtional packages can be installed in spyder-env e.g. with conda 

**Installing setup.py for Windows
- opening setup.oy via IDLE (Python GUI)
- starting the program via Run ... Customized (Shift + F5)
- run program with the command "install"
