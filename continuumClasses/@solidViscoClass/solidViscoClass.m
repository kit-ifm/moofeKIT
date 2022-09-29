classdef solidViscoClass < solidSuperClass  
    properties
        epsilonViscoN
        epsilonViscoN1
        dissWorkN
        dissWorkN1
        alpha                                                              
        r                                                                   
        p0
        pInf
        eRef                                                             
        eps0
        sigma
        epsilon
        totalTimeSteps
    end
    properties (Constant = true)
        additionalFields = 0;
    end
    methods
        %% constructor
        function obj = solidViscoClass(dofObject)
            if nargin==0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
        end
    end
    
    methods
        %% mandatory methods
        function updateHistoryField(obj,dofObject)
            obj.epsilonViscoN = obj.epsilonViscoN1;
        end
        function initializeEpsilonViscoN1(obj)
            if isempty(obj.epsilonViscoN1)
                obj.alpha=0.4;
                obj.p0= (1-obj.alpha)*obj.materialObject.eModul0;
                obj.r=0.5;
                obj.eps0=0.4;
                obj.pInf=(obj.alpha)*obj.materialObject.eModul0;
                obj.eRef=obj.materialObject.eModul1; 
                numberOfElements = size(obj.meshObject.edof,1);
                numberOfGausspoints = obj.shapeFunctionObject.numberOfGausspoints;
                obj.epsilonViscoN1 = zeros(numberOfElements, numberOfGausspoints);
            end
            if isempty(obj.dissWorkN)
                obj.dissWorkN=0;
            end            
            if isempty(obj.dissWorkN1)
                obj.dissWorkN1=0;
            end
            if isempty(obj.sigma)
                 obj.sigma=zeros(obj.totalTimeSteps,1);
                 obj.epsilon=zeros(obj.totalTimeSteps,1);  
            end
        end
        
 
    end
end