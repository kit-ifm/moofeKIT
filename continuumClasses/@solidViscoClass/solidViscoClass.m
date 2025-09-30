classdef solidViscoClass < solidSuperClass
    properties
        flagHistoryFields = true;
        epsilonViscoN % internal viscoelastic variable at timestep N
        epsilonViscoN1 % internal viscoelastic variable at timestep N+1
        dissWorkN % sum of dissipated work over time t_0 to t_n
        dissWorkN1 % sum of dissipated work over time t_0 to t_n+1
        linearity = 'linear' % switch for linear or nonlinear material behaviour
        alpha % ?
        r % ?
        p0 % ?
        pInf % ?
        eRef % ?
        eps0 % ?
        sigma % should not be property of class?
        epsilon % should not be property of class?
        totalTimeSteps % already defined in setupObject?
    end
    properties
        additionalFields = 0;
        dofsPerAdditionalField = [];
    end
    methods

        %% constructor
        function obj = solidViscoClass(dofObject)
            if nargin == 0
                error('dofObject is required as input');
            end
            dofObject.listContinuumObjects{end+1} = obj;
        end
    end

    methods

        %% mandatory methods
        function set.linearity(obj, value)
            if strcmpi(value, 'linear') || strcmpi(value, 'nonlinear')
                obj.linearity = value;
            else
                error('linearity property must be linear or nonlinear')
            end
        end
        function updateHistoryField(obj, dofObject)
            obj.epsilonViscoN = obj.epsilonViscoN1;
        end
        function initializeEpsilonViscoN1(obj)
            if isempty(obj.epsilonViscoN1)
                numberOfElements = size(obj.meshObject.edof, 1);
                numberOfGausspoints = obj.shapeFunctionObject.numberOfGausspoints;
                dimension = obj.dimension;
                if strcmp(obj.linearity, 'linear')
                    % linear case: scalar internal variabel on every Gaussian point
                    obj.alpha = 0.4;
                    obj.p0 = (1 - obj.alpha) * obj.materialObject.eModul0;
                    obj.r = 0.5;
                    obj.eps0 = 0.4;
                    obj.pInf = (obj.alpha) * obj.materialObject.eModul0;
                    obj.eRef = obj.materialObject.eModul1;
                    obj.epsilonViscoN1 = zeros(numberOfElements, numberOfGausspoints);
                else %strcmp(obj.linearity,'nonlinear')
                    % nonlinear case: internal variabel: symmetric tensor C_i (Voigt Notation) on every element
                    obj.epsilonViscoN1 = zeros(numberOfElements, numberOfGausspoints, sum(1:dimension));
                    for e = 1:numberOfElements
                        for k = 1:numberOfGausspoints
                            obj.epsilonViscoN1(e,k,:)=[1, 1, 1, 0, 0, 0];
                        end
                    end
                end
            end
            if isempty(obj.dissWorkN)
                obj.dissWorkN = 0;
            end
            if isempty(obj.dissWorkN1)
                obj.dissWorkN1 = 0;
            end
            if isempty(obj.sigma)
                obj.sigma = zeros(obj.totalTimeSteps, sum(1:dimension)); % see properties!
            end
            if isempty(obj.epsilon)
                obj.epsilon = zeros(obj.totalTimeSteps, sum(1:dimension)); % see properties!
            end
        end


    end
end