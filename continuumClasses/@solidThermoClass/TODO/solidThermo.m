classdef solidThermo < solver.conti
    %% Class definition for thermo solids (Lagrange, NURBS, HNURBS)
    %
    % Created: 			Mi, 10 Aug 2016
    % Responsibilty: 	Mark Schiebl
    % Editors:
    % g005540  gk538  gk604  gk612
    % Description:
    %----------------------------------------------------------------------
    % TODO Mark:    1. GENERIC Routines to FE-Routines
    %               2. ProjectThermo Routines to FE-Routines
    %
    
    properties
        PLOTData            % Data for postprocessing
        ProjVsave           % Projection Values
        QREF                % Nodel geometry matrix - reference configuration.
        QN                  % Nodel geometry matrix at time n.
        QN1                 % Nodel geometry matrix at time n+1.
        VN                  % Velocity  at time n.
        PN                  % Linear mometum at time n.
        MAT                 % Material data.
        BC                  % Boundary Conditions (Zero important for staggered routine)
        EPOT = 0;           % Strain energy of the system.
        helmholtz = 0;      % free helmholtz energy of the system
        e = 0;              % internal energy of the system
        entr = 0;           % entropy of the system
        GRAVITY             % Gravity vector.
        ROUTINE             % Element routine to be evaluated.
        REFTEMP             % Reference temperatur.
        WEIGHTS             % Weights of Nodes (NURBS)
        ENH = false;        % Logical, true if enhanced modes are used.
        ENHStruct           % Structure containing information to recover enhanced modes.
        ALPHAN              % Enhanced modes at time n.
        ALPHAN1             % Enhanced modes at time n+1.
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
        TAU                 % Generic in thermovariable tau in [theta,entropy,energy,phi(tau)]
        INITTAU;            % Generic initial values in thermovariable inittau in [theta,entropy,energy,phi(tau)]
        INITTRANS = false;  % true if transformation of initial condition in generic format is required
        numtang = false;    % Numerical Tangenent
        L2Proj              % GENERIC in Entropy, L2 Proj = True / False
        INITTRANSPATH       % Transformationpath for GENERIC
        ORDER_MIXED_VAR
        SHAPEF_LM           % Shape functions for LM field
        StaConStruct        % Structure containing information for static condensation process.
        SIGMAN          % Enhanced stresses at time n.
        SIGMAN1         % Enhanced stresses at time n.
        EPSILONN        % Enhanced strains at time n.
        EPSILONN1       % Enhanced strain at time n.
        SIGMA_FN        % Work Conjugated Stress (w.r.t Deformation Gradient F) at time n
        SIGMA_FN1       % Work Conjugated Stress (w.r.t Deformation Gradient F) at time n1
        SIGMA_HN        % Work Conjugated Stress (w.r.t Area map tensor H) at time n
        SIGMA_HN1       % Work Conjugated Stress (w.r.t Area map tensor H) at time n1
        SIGMA_JN        % Work Conjugated Stress (w.r.t Determinant) at time n
        SIGMA_JN1       % Work Conjugated Stress (w.r.t Determinant) at time n1
        FN              % Independent Deformation Gradient F at time n
        FN1             % Independent Deformation Gradient F at time n1
        HN              % Independent Area map tensor H at time n
        HN1             % Independent Area map tensor H at time n1
        JNe             % Independent Jacobian at time n
        JN1e            % Independent Jacobian at time n1
        SHAPEF_F        % Shape functions for F (mixed methods)
        SHAPEF_Cof      % Shape functions for Cofactor (mixed methods)
        SHAPEF_Det      % Shape functions for Determinant (mixed methods)
        ModData         % Modification Data for basis modification
        ModificatedElements
    end
    properties (Dependent = true)
        JN           % Total angular momentum at time n.
        L2ProjV      % Node Values of L2Proj
    end
    
    methods
        %% Constructor
        function obj = solidThermo(mesh,varargin)
            if nargin == 0
            else
                obj = constructor(obj,mesh,varargin{:});
                if numel(obj) == 1 && isempty(obj(1).ROUTINE)
                    obj(1) = [];
                end
            end
        end
        
        %% FE- Routines
        %FIXME: Nicht alle Routinen aufgefuehrt
        out = disp_ogdenEndPoint(obj,varargin)                   % Ogden strain-energy function, evaluated at the end-point of the configuration.
        out = disp_ogdenDiscreteGradient(obj,varargin)           % Ogden strain-energy function, using the discrete gradient.
        out = disp_ogdenMixed(obj,varargin)                      % Ogden strain-energy function, evaluated at the mid-point of the strains.
        out = disp_ogdenMidPoint(obj,varargin)                   % Ogden strain-energy function, evaluated at the mid-point of the configuration.
        out = disp_ogdenEnhDiscreteGradient(obj,varargin)        % Ogden strain-energy function, using the discrete gradient and enhanced modes.
        out = disp_ogdenEnhMixed(obj,varargin)                   % Ogden strain-energy function, evaluated at the mid-point of the strains and using enhanced modes.
        out = disp_ogdenEnhMidPoint(obj,varargin)                % Ogden strain-energy function, evaluated at the mid-point of the configuration and using enhanced modes.
        out = disp_NeoHookEndPoint(obj,varargin)                 % NeoHook strain-energy function, evaluated at the end-point of the configuration.
        out = disp_NeoHookMixed(obj,varargin)                    % NeoHook strain-energy function, evaluated at the mid-point of the strains.
        out = disp_NeoHookMidPoint(obj,varargin)                 % NeoHook strain-energy function, evaluated at the mid-point of the configuration.
        out = disp_MooneyRivlinEndPoint(obj,varargin)            % MooneyRivlin strain-energy function, evaluated at the end-point of the configuration.
        out = wedge_MooneyRivlinEndPoint(obj,varargin)           % MooneyRivlin wedge strain-energy function, evaluated at the end-point of the configuration.
        out = mmhr_MooneyRivlinEndPoint(obj,varargin)            % MooneyRivlin hellinger-reissner functional, evaluated at the end-point of the configuration.
        out = generic_neoHooke_entropy_midPoint(obj,varargin)    % Generic in Entropy form for NeoHook strain-energy function, evaluated at the mid-point configuration.
        out = generic_neoHooke_entropy_L2_midPoint(obj,varargin) % Generic in Entropy form for NeoHook strain-energy function with L2 Projection field, evaluated at the mid-point configuration.
        out = generic_neoHooke_temp_midPoint(obj,varargin)       % Generic in Temperature form for NeoHook strain-energy function, evaluated at the mid-point configuration.
        out = generic_neoHooke_temp_midPoint_L2(obj,varargin)    % Generic in Temperature form for NeoHook strain-energy function with L2 Projection field, evaluated at the mid-point configuration.
        out = generic_neoHooke_entropy_TC(obj,varargin)          % Generic in Entropy form for NeoHook strain-energy function, TC integrator (not true TC, just Entropy consistent)-
        out = generic_neoHooke_entropy_TC_L2(obj,varargin)       % Generic in Entropy form for NeoHook strain-energy function with L2 Projection, TC integrator-
        out = generic_neoHooke_temp_TC(obj,varargin)             % Generic in Temperature form for NeoHook strain-energy function, TC integrator (not true TC, just Energy consistent).
        out = generic_neoHooke_temp_TC_L2(obj,varargin)          % Generic in Temperature form for NeoHook strain-energy function, TC integrator.
        out = generic_wedge_neoHooke_temp_TC_L2(obj,varargin)    % Generic in Temperature form for NeoHook strain-energy function using Tensor cross notion, TC integrator.
        
        %% Tools
        obj = constructor(obj,mesh)                         % Constructor.
        obj = selectRoutine(obj,integrator)                 % Select element routine.
        obj = update(obj,varargin)                          % Updates specified field
        plot(obj,varargin)                                  % Plot tool.
        
        %% Get and set functions
        function out = get.BC(obj)
            %dirichlet BC for staggered mode
            
            DIM = obj.DIM;
            ACTIVEFIELD = obj.ACTIVEFIELD;
            if ACTIVEFIELD == 1
                list = DIM+1:size(obj.NODESDOF,2);
            elseif ACTIVEFIELD > 1
                list = 1:size(obj.NODESDOF,2);
                list(DIM+ACTIVEFIELD-1) = [];
            else
                list = [];
            end
            zw = obj.NODESDOF(:,list);
            out = unique(zw(:));
        end
        
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
        function out = get.L2ProjV(obj)
            if obj.L2Proj == true
                fdof = obj.NODESDOF(:,obj.DIM+1);
                maxElAct = max(fdof);
                internal=residual(obj,'L2Proj');
                rHam = full(sparse(vertcat(internal.structure(:).edofE),1,vertcat(internal.structure(:).Re),maxElAct,1));
                kHam = sparse(vertcat(internal.structure(:).pI), vertcat(internal.structure(:).pJ), vertcat(internal.structure(:).pK),maxElAct,maxElAct);
                out = kHam(fdof,fdof)\rHam(fdof,1);
            else
                error('No L2 Projection necessary')
            end
        end
        function out = get.PLOTData(obj)
            if strcmp(obj.DESIGNATOR,'generic')
                if strcmp(obj.TAU,'entr')
                    fdof = obj.NODESDOF(:,obj.DIM+obj.ADDFIELDS);
                    maxElAct = max(fdof);
                    internal=residual(obj,'L2Proj');
                    rHam = full(sparse(vertcat(internal.structure(:).edofE),1,vertcat(internal.structure(:).Re),maxElAct,1));
                    kHam = sparse(vertcat(internal.structure(:).pI), vertcat(internal.structure(:).pJ), vertcat(internal.structure(:).pK),maxElAct,maxElAct);
                    out = kHam(fdof,fdof)\rHam(fdof,1);
                elseif strcmp(obj.TAU,'engy')
                    fdof = obj.NODESDOF(:,obj.DIM+obj.ADDFIELDS);
                    maxElAct = max(fdof);
                    internal=residual(obj,'L2Proj');
                    rHam = full(sparse(vertcat(internal.structure(:).edofE),1,vertcat(internal.structure(:).Re),maxElAct,1));
                    kHam = sparse(vertcat(internal.structure(:).pI), vertcat(internal.structure(:).pJ), vertcat(internal.structure(:).pK),maxElAct,maxElAct);
                    out = kHam(fdof,fdof)\rHam(fdof,1);
                end
            else
                error('Ploterror');
            end
        end
    end
end
