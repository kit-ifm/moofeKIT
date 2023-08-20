# moofeKIT (matlab object-oriented finite element kit)

### changelog
https://git.scc.kit.edu/ifm/moofeKIT/-/blob/master/changelog.md
- meetings are scheduled approximately every three weeks (the upcoming meeting is scheduled in the changelog)
- the users are encouraged to merge all feature branches (**if well-engineered, useful to others and intended for the public**) to masterExperimental branch generously before the upcoming meeting 
- masterExperimental branch will be merged to master branch
- main changes will be logged
- **only the master branch will be merged to master branch of public github repo** https://github.com/kit-ifm/moofeKIT
### git howto and workflow idea
git howto and workflow idea is provided via https://git.scc.kit.edu/ifm/anleitungen/software/-/blob/master/git/simple_cheat_sheet.md

see also https://training.github.com/downloads/github-git-cheat-sheet.pdf and https://rhoenerlebnis.de/_upl/de/_pdf-seite/git_cheatsheet_de_white.pdf
### best practice for clean code
- Language should be english for everything (variables, classes, folders, files, comments etc.).
- Write complete names for variables, classes, folders, files etc. with convention: 
  - First letter of first name should be written in lower case. First letter of (second, third etc) joined words should be capitalized.
  - Example >> parameterMaterialNu = 3; (not like >> paramnu = 3;)
- Do not duplicate code (minimalism)!
- A function should do one thing, and only one thing (modularity)!
- Use external toolbox https://github.com/davidvarga/MBeautifier for source code formatting. It is already included via submodule 'git submodule add https://github.com/davidvarga/MBeautifier.git externalToolboxes/MBeautifier'. In matlab the code in current editor can be 'beautified' via command 'MBeautify.formatCurrentEditorPage()'. A shortcut can be created via 'Home', 'Favorites', 'New Favorite' by adding the command 'MBeautify.formatCurrentEditorPage()'.
### handling issues
- please list issues https://git.scc.kit.edu/franke/moofeKIT/-/issues
- if the problem could be devided into small parts **use checklist** in markdown 
- assign it (if possible) to a member
- set due date (if possible)
### howto comment/manual
a dedicated manual is not planned, all scritps/classes/functions/methods should be documented via commented help text immediately below the function definition line, see 
https://de.mathworks.com/help/matlab/matlab_prog/add-help-for-your-program.html
- example files
  - script: https://git.scc.kit.edu/ifm/moofeKIT/-/blob/master/scripts/pre/LShapeElectroThermoMechanics.m
  - class: https://git.scc.kit.edu/ifm/moofeKIT/-/blob/master/coreClasses/@solidSuperClass/solidSuperClass.m
  - function: https://git.scc.kit.edu/ifm/moofeKIT/-/blob/master/continuumClasses/@solidClass/displacementSCSaintVenantEndpoint.m
- useful matlab code
  - for help of a function/class type _>> help function-name_
                                   or _>> doc function-name_
  - for look for a keyword _>> lookfor keyword_
  - location of filename _>> which filename_ 
  - list properties of a class name _>> properties className_
  - list methods of a class name _>> methods className_
  - for TODO/FIXME report go to current folder browser and then select Reports > TODO/FIXME Report
  
some additional manual-sheets can be found in folder manual/
### object-oriented programming
https://en.wikipedia.org/wiki/Object-oriented_programming or german version https://de.wikipedia.org/wiki/Objektorientierte_Programmierung

matlab specific object-oriented programming can be found at https://de.mathworks.com/products/matlab/object-oriented-programming.html

for a simple matlab example https://de.mathworks.com/company/newsletters/articles/introduction-to-object-oriented-programming-in-matlab.html might be helpful
### gitlab runner for continuous integration
a short manual for gitlab runner for continous integration (CI) with testing code in matlab can be found at https://git.scc.kit.edu/franke/ubuntuinstallationmitarbeiterrechnerifm/-/blob/master/gitlabRunner.txt
