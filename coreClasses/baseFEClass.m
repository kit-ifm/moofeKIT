classdef (Abstract) baseFEClass < matlab.mixin.Copyable
% BASEFECLASS Parent class for all finite element subclasses which defines 
% necessary methods for all subclasses.
% 
% INHERIT
% Inherit a child class via classdef childClass < baseFEClass. Note, an
% abstract class can not be initiated. 
% 
% GENERAL
% An abstract class serves as a basis (that is, a superclass) for a group 
% of related subclasses. An abstract class can define abstract properties 
% and methods that subclasses implement. Each subclass can implement the 
% concrete properties and methods in a way that supports their specific re-
% quirements.
% Instead of a handle class a copyable class is taken:
% matlab.mixin.Copyable look at https://de.mathworks.com/help/matlab/matlab_oop/custom-copy-behavior.html
% for further informations.
% 
% REFERENCE(S)
% https://de.mathworks.com/help/matlab/matlab_oop/abstract-classes-and-interfaces.html
% 
% SEE ALSO
% solidSuperClass, 
% dirichletClass, 
% dirichletLagrangeClass, 
% neumannClass
% 
% CREATOR(S) 
% Marlon Franke

    properties (Abstract = true, Constant = true)
        callMassMatrix
        callElements
    end
    properties
        meshObject
    end
    methods
        function obj = baseFEClass()
            obj.meshObject = meshClass();
        end
    end
    methods (Abstract = true)
        initializeShapeFunctions(obj)
        initializeGlobalDofs(obj)
        initializeQV(obj)
        updateGlobalField(obj)
        updateContinuumFieldPreNewtonLoop(obj)
        updateTimeDependentFieldPreNewtonLoop(obj)
        updateTimeDependentFieldNewtonLoop(obj)
        updateContinuumFieldNewtonLoop(obj)
        updateContinuumFieldPostNewtonLoop(obj)
    end
end