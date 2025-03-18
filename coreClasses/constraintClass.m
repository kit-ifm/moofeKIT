classdef constraintClass < handle
    %CONSTRAINTCLASS 
    
    properties
        dofObject
        constrainedElement
        constrainedDOFinElement
        constrainedDOFinOtherElement
    end
    
    methods
        %%%%%%%%%%%%%%%%%
        function obj = constraintClass(dofObject)
            obj.dofObject = dofObject;
            dofObject.isConstrained = true;
            dofObject.constraintObjectList{end+1} = obj;
        end

        function [dataFE, doSolve] = implementConstraints(obj, dataFE, doSolve)
                
            for constrainedCounter = 1:numel(obj.constrainedDOFinElement)
                constrainedDOF = obj.constrainedDOFinElement(constrainedCounter);
                indexReI = dataFE(obj.constrainedElement).indexReI==constrainedDOF;
                indexKeI = dataFE(obj.constrainedElement).indexKeI==constrainedDOF;
                indexKeJ = dataFE(obj.constrainedElement).indexKeJ==constrainedDOF;


                dataFE(obj.constrainedElement).indexReI(indexReI) = obj.constrainedDOFinOtherElement(constrainedCounter);
                dataFE(obj.constrainedElement).indexKeI(indexKeI) = obj.constrainedDOFinOtherElement(constrainedCounter);
                dataFE(obj.constrainedElement).indexKeJ(indexKeJ) = obj.constrainedDOFinOtherElement(constrainedCounter);

                doSolve(obj.constrainedDOFinElement) = 0;
            end

        end

        function massDataFE = implementConstraintsInMass(obj, massDataFE)
                
            for constrainedCounter = 1:numel(obj.constrainedDOFinElement)
                constrainedDOF = obj.constrainedDOFinElement(constrainedCounter);
                indexMeI = massDataFE(obj.constrainedElement).indexMeI==constrainedDOF;
                indexMeJ = massDataFE(obj.constrainedElement).indexMeJ==constrainedDOF;


                massDataFE(obj.constrainedElement).indexMeI(indexMeI) = obj.constrainedDOFinOtherElement(constrainedCounter);
                massDataFE(obj.constrainedElement).indexMeJ(indexMeJ) = obj.constrainedDOFinOtherElement(constrainedCounter);

            end

        end

        function updateConstrainedDOF(obj)
            obj.dofObject.qN1(obj.constrainedDOFinElement) = obj.dofObject.qN1(obj.constrainedDOFinOtherElement);
        end
    end
end

