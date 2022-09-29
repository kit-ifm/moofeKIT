classdef solidElectroThermo < solver.conti
    %% Class definition for electrothermo solids (Lagrange, NURBS or HNURBS).
    %
    % 08.04.2019 M. Franke
    
    properties
        QREF                % Nodel geometry matrix - reference configuration.
        QN                  % Nodel geometry matrix at time n.
        QN1                 % Nodel geometry matrix at time n+1.
        VN                  % Velocity  at time n.
        PN                  % Linear mometum at time n.
        MAT                 % Material data.
        BC = [];            % Boundary Conditions (Zero important for staggered routine)
        EPOT = 0;           % Strain energy of the system.
        helmholtz = 0;      % free helmholtz energy of the system
        e = 0;              % internal energy of the system
        entr = 0;           % entropy of the system
        GRAVITY             % Gravity vector.
        ROUTINE             % Element routine to be evaluated.
        REFPHI              % Reference potential.
        REFTEMP             % Reference temperatur.
        WEIGHTS             % Weights of Nodes (NURBS)
        HCONTROL            % Controls h-refinement.
        levelVec
        levelStruct
        IDSTRUCT
        IDTOPOLOGIE
        STMatrix
        ELEMENTS
        KNOTS
        ORDER
        DESIGNATOR          % Identifier for routines.
        FN                  % Independent Deformation Gradient F at time n
        FN1                 % Independent Deformation Gradient F at time n1
        HN                  % Independent Area map tensor H at time n
        HN1                 % Independent Area map tensor H at time n1
        JNe                 % Independent Jacobian at time n
        JN1e                % Independent Jacobian at time n1  
        DN                  % Independent electrical displacement (material) at tn
        DN1                 % Independent electrical displacement (material) at tn1
        dN                  % Independent electrical displacement (spatial) at tn
        dN1                 % Independent electrical displacement (spatial) at tn1
        SIGMA_FN            % Work Conjugated Stress (w.r.t Deformation Gradient F) at time n
        SIGMA_FN1           % Work Conjugated Stress (w.r.t Deformation Gradient F) at time n1
        SIGMA_HN            % Work Conjugated Stress (w.r.t Area map tensor H) at time n
        SIGMA_HN1           % Work Conjugated Stress (w.r.t Area map tensor H) at time n1
        SIGMA_JN            % Work Conjugated Stress (w.r.t Determinant) at time n
        SIGMA_JN1           % Work Conjugated Stress (w.r.t Determinant) at time n1
        SIGMA_DN            % Work Conjugated Stress w.r.t electrical displacement (material) at tn
        SIGMA_DN1           % Work Conjugated Stress w.r.t electrical displacement (material) at tn1
        SIGMA_dN            % Work Conjugated Stress w.r.t electrical displacement (spatial) at tn
        SIGMA_dN1           % Work Conjugated Stress w.r.t electrical displacement (spatial) at tn1
        SHAPEF_F            % Shape functions for F 
        SHAPEF_Cof          % Shape functions for Cofactor 
        SHAPEF_Det          % Shape functions for Determinant
        SHAPEF_D            % Shape functions for electrical displacement (material) 
        SHAPEF_d            % Shape functions for electrical displacement (spatial) 
        SHAPEF_phi          % Shape functions for electrical potential 
        ORDER_MIXED_VAR     % Interpolation order of F, Cof, Det
        ORDER_MIXED_VAR_EL  % Interpolation order of D ,d  (electrical displacement)
        StaConStruct        % Structure containing information for static condensation process.
        numericalTangent = false;
    end
    properties (Dependent = true)
        JN           % Total angular momentum at time n.
    end
    
    methods
        %% Constructor
        function obj = solidElectroThermo(mesh,varargin)
            if nargin == 0
            else
                obj = constructor(obj,mesh,varargin{:});
                if numel(obj) == 1 && isempty(obj(1).ROUTINE)
                    obj(1) = [];
                end
            end
        end
        
        %% FE- Routines
        out = mmhwSC_mooneyRivlin_endPoint(obj,varargin)           % MooneyRivlin hu-washizu functional, evaluated at the end-point of the configuration.

        %% Tools
        obj = constructor(obj,mesh)                         % Constructor.
        obj = selectRoutine(obj,integrator)                 % Select element routine.
        obj = update(obj,varargin)                          % Updates specified field
        plot(obj,varargin)                                  % Plot tool.
        function out = get.JN(obj)                          % Get method for the angular momentum.
            out = 0;
            if size(obj.QN,2) == 0 || size(obj.PN,2) == 0
                if obj.DIM == 2
                    out = 0;
                elseif obj.DIM == 3
                    out = [0 0 0]';
                end
            else
                if obj.DIM == 2
                    out = sum(obj.QN(:,1).*obj.PN(:,2) - obj.QN(:,2).*obj.PN(:,1));
                elseif obj.DIM == 3
                    out = [sum(obj.QN(:,2).*obj.PN(:,3) - obj.QN(:,3).*obj.PN(:,2));
                        sum(obj.QN(:,3).*obj.PN(:,1) - obj.QN(:,1).*obj.PN(:,3));
                        sum(obj.QN(:,1).*obj.PN(:,2) - obj.QN(:,2).*obj.PN(:,1))];
                end
            end
        end
end
end
