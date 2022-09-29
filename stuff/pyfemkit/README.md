# finiteElementCode

**Install module/package pyfemkit**  
pip install -e .

**TODO later (as soon as framework is fixed)**
- [ ] define testroutine implement git-CI (py and matlab-runner)
- [ ] basic feature: universal numerical tangent (FD & complex tangent)
- [ ] mixed finite elements -> solver/static condensation
- [ ] implement GFE routines
- [ ] implement FEidF routines
- [ ] implement KEuG routines
- [ ] implement abaqus-input-file interface

**The Zen of Python, by Tim Peters**

Beautiful is better than ugly.  
Explicit is better than implicit.  
Simple is better than complex.  
Complex is better than complicated.  
Flat is better than nested.  
Sparse is better than dense.  
Readability counts.  
Special cases aren't special enough to break the rules.  
Although practicality beats purity.  
Errors should never pass silently.  
Unless explicitly silenced.  
In the face of ambiguity, refuse the temptation to guess.  
There should be one-- and preferably only one --obvious way to do it.  
Although that way may not be obvious at first unless you're Dutch.  
Now is better than never.  
Although never is often better than *right* now.  
If the implementation is hard to explain, it's a bad idea.  
If the implementation is easy to explain, it may be a good idea.  
Namespaces are one honking great idea -- let's do more of those!  

**Best practice for clean code**
- Language should be english for everything (variables, classes, folders, files, comments etc.).
- Write complete names for variables, classes, folders, files etc. with convention: 
  - First letter of first name should be written in lower case. First letter of (second, third etc) joined words should be capitalized.
  - Example >> parameterMaterialNu = 3; (not like >> paramnu = 3;)
- Do not duplicate code (minimalism)!
- A function should do one thing, and only one thing (modularity)!
