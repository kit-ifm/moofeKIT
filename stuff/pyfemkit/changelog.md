# changelog
**30/04/2021**
- Python
    - [x] object-oriented FE code
    - [x] paralellization (mass matrix not yet)
        - 'Singlecore' Python ~1.5x slower (Matlab with paralellization of direct solver)
        - 'Multicore' (24 cores): 13s (Matlab) to 16s (Python without paralellization of mass matrix)
    - [ ] todo?: automatic differentiation

**19/03/2021**
- Matlab
    - [x] core change: object search abolished, catch objects via property of classes by creation of an instance implemented via idea of Robin
    - [x] paraview (and blender) interface (from Rogelio)
        - VTKplot-routine
        - automatic call of paraview
        - png/eps etc. export
        - generate movies
- Python
    - [x] diagrams (2D)
        - plot with matplotlib
        - tikz export via tikzplotlib

**12/03/2021**
- Matlab
    - core changes
        - [x] auto updates via dofObject
        - [x] flags in classes and runNewton-script for update-procedures (idea: for further classes hopefully no need for adjustement of runNewton-script)
        - [ ] failed: catch objects via property of classes by creation of an instance
        - [x] baseFEClass (nodes, edof, globalFullEdof, globalNodesDof, dimension)
        - [x] solidSuperClass (qR, qN1, qN, vN1, vN) (derived from baseFEClass)
        - [x] clear/destruct all 'pre' Objects except continuumObjects
        - [x] function structure of runNewton-script
        - [x] speed up: save dofObject -> object search
    - proof of concept
        - [x] discrete Gradient Elementroutine
        - [x] solidThermoClass (proof of concept for auto updates)
        - [x] numerical Tangent (FD) function-based structure (parfor-loop not possible)
    - additional changes
        - [x] dofs -> dof (e.g. dofsObject -> dofObject)
        - [x] dirichletClass -> implement full features (arbritary dof, timeFunction: include nodal information)

**02/03/2021**
- Matlab (first implementation of the code, idea is to determine the framework in Matlab and afterwards to test a Python implementation)
    - structure of the folders (coreClasses, scripts, commonFunctions, meshes, continuumClasses, etcs)
    - code (basic elements/classes)
        - [x] prepareWorkspace
            - [x] setupClass
            - [x] solidClass
                - edof -> local Element-Knotenfreiheitsgradtabelle
                - (fullEdof -> local, full element-degree-of-freedom table)
                - (globalEdof -> global element-degree-of-freedom table, needed for more than one continuumObject)
                - globalFullEdof -> global, full element-degree-of-freedom table
                - (nodesDof -> local full nodal degrees-of-freedom list)
                - globalNodesDof -> global full element-degree-of-freedom table
            - [x] dirichletClass
            - [x] runNewton
- Python 
    - Spyder-IDE
        - editor
        - command window -> console
        - variable explorer
        - plots
        - profiler and code analyzer
    - [x] python implementation of the code of basic fe course ('Lehrcode')
        - [x] suplementation of sparse assembly and paralellization
        - [x] symbolic computing
