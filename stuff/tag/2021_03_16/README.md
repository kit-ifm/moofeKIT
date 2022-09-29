# finiteElementCode

**changelog 12/03/2021**
- core changes
    - [x] auto updates
    - [x] flags in classes and runNewton-script for update-procedures (idea: for further classes hopefully no need for adjustement of runNewton-script)
    - [failed] catch objects via property of classes by creation of an instance
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

**TODO**
- [] basic feature: universal numerical tangent (FD & complex tangent)
- [] proof of concept: mixed finite elements -> solver/static condensation
- [] test python implementation of the above

**TODO later (as soon as framework is fixed)**
- [] define testroutine implement git-CI (matlab-runner)
- [] implement GFE routines
- [] implement FEidF routines
- [] implement KEuG routines
- [] implement abaqus-input-file interface
- [] implement class for startScript config
- [] evaluate automatic classtree for all classes (including inherited abstract and/or superclasses)

**Best practice for clean code**
- Language should be english for everything (variables, classes, folders, files, comments etc.).
- Write complete names for variables, classes, folders, files etc. with convention: 
  - First letter of first name should be written in lower case. First letter of (second, third etc) joined words should be capitalized.
  - Example >> parameterMaterialNu = 3; (not like >> paramnu = 3;)
- Do not duplicate code (minimalism)!
- A function should do one thing, and only one thing (modularity)!
