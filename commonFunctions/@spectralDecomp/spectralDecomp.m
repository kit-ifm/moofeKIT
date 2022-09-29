classdef spectralDecomp < handle
    %SPECTRALDECOMP 
    %class to compute a spectral decomposition 
    %Robin Pfefferkorn 30.8.2019
    
    %eigenvalue based formulation:
    %----------------------------------------------------------
    %-various algorithms 
    %-spatial as well as material formulation
    
    %changelog
    %----------------------------------------------------------
    %v1.0 (30.8.2019) 
    %first release
    %
    %v1.1 (21.10.2020) 
    %added version without perturbation using exact eigen vectors and
    %adjusted tangent expressions in case of identical eigenvalues
    
    %##################################################################
    %Properties
    %##################################################################
    properties
        %general settings concerning the algorithm
        algorithm = 'perturbation';
        perturbThreshold = 1e5*eps; %difference of eigenvalues at which they are considered equal
        %recommended values:
        %perturbation algorithm: 1e5*eps
        %exactEigenVec algorithm: 1e2*sqrt(eps) (see Jeremic, Cheng (2005))
        spatialFormulation = true;
    end
    properties(SetAccess=private)
        %general settings concerning the algorithm
        algorithmNr = 0;
        %eigenvalues and bases
        principalStretch = zeros(3,1) %column vector containing the 3 principal stretches
        eigenVectors = zeros(3,3); %matrix containing all 3 eigenvectors in the columns
        eigenBases = cell(6,1); %cell containing the 6 eigenvector bases
        eigenBasesVoigt = []; %2D matrix containing the 6 eigenvector bases in voigt notation
    end
    
    %##################################################################
    %public methods
    %##################################################################
    methods
        %------------------------------------------------------------------
        %constructor
        %------------------------------------------------------------------
        function obj = spectralDecomp(algorithm,spatial)
            if nargin==0
                %do nothing (use std values)
            elseif nargin==1
                obj.algorithm = algorithm;
            elseif nargin==2
                obj.algorithm = algorithm;
                obj.spatialFormulation = spatial;
            else
                error('wrong amount of input arguments')
            end
        end
        %------------------------------------------------------------------
        %properties get set
        %------------------------------------------------------------------
        function set.algorithm(obj,inp)
            if ~any(strcmpi(inp,{'perturbation','exacteigenvec'}))
                error('algorithm may only be: perturbation, exacteigenvec')
            else
                obj.algorithm = inp;
                setAlgorithmNr(obj);
            end
        end
        function set.perturbThreshold(obj,inp)
            if inp<=0
                error('perturbThreshold must be >0')
            else
                obj.perturbThreshold = inp;
            end
        end
        function set.spatialFormulation(obj,inp)
            if ~islogical(inp)
                error('spatialFormulation must be true/false')
            else
                obj.spatialFormulation = inp;
            end
        end
        
        %------------------------------------------------------------------
        %further functions (actual computations)
        %------------------------------------------------------------------
        obj=decomposition(obj,b,mapVoigt)
        out=assemTens(obj,lam)
        out=assemMatTang(obj,lam,eMat,factStress)
        obj=scaleEigenvectors(obj,factors)
    end
    
    %##################################################################
    %private methods
    %##################################################################
    methods(Access=private,Hidden=true)
        function setAlgorithmNr(obj)
            %sets the algorithm number (faster comparison)
            switch lower(obj.algorithm)
                case 'perturbation'
                    obj.algorithmNr = 0;
                case 'exacteigenvec'
                    obj.algorithmNr = 1;
                otherwise
                    error('not implemented')
            end
        end
    end
end

