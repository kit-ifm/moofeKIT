function displacementHyperelasticMidpoint(obj,setupObject,varargin)
%% Creates the residual and the tangent of the given obj.
%
% out = disp_logStrain_endPoint(obj,'PropertyName',PropertyValue)
%
% Description:
% -various hyperelastic laws,
% -evaluated at time n+1/2, i.e. midpoint rule.
% -spatial formulation in Tau and left Cauchy Green 
%
% 19.4.2021 Robin Pfefferkorn

%% Check input
computeStresses = 0;%stresses for postprocessing. if computeStresses=0 -> compute standard residual and tangent
ii = 1;
while ii <= size(varargin,2)
    if strcmpi(varargin{ii},'computeStresses')
        computeStresses = varargin{ii+1};
        ii = ii+1;
    end
    ii = ii+1;
end
%% setup
% load objects
meshObject = obj.meshObject;
shapeFunctionObject = obj.shapeFunctionObject;
storageFEObject = obj.storageFEObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
%   mixedFEObject = obj.mixedFEObject;
% element degree of freedom tables and more
edof = meshObject.edof;
globalFullEdof = meshObject.globalFullEdof;
numberOfElements = size(globalFullEdof,1);
numberOfDofs = size(globalFullEdof,2);
dimension = obj.dimension;
% gauss integration and shape function
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
% nodal dofs
qR = obj.qR;
qN1 = obj.qN1;
qN = obj.qN;
qN05 = 0.5*(qN + qN1);
% material data and voigt notation
I = eye(3);
selectMapVoigt(mapVoigtObject,dimension,'symmetric');
% initialization energy and storageFEObject
strainEnergy = 0;
initializeDataFE(storageFEObject);
%% element loop for residual and tangent
for e = 1:numberOfElements
    %for eigenvalue materials. test if stresses are realistic
    stretchOutOfRange = false;
    % Element routine (residual and tangent)
    Re = zeros(numberOfDofs,1);
    Ke = zeros(numberOfDofs);
    %post processing stresses
    se = zeros(numberOfDofs/dimension,1);
    Me = zeros(numberOfDofs/dimension);
        
    %dofs and jacobian
    edN1 = qN1(edof(e,:),1:dimension).'; %#ok<PFBNS>
    edN05 = qN05(edof(e,:),1:dimension).'; %#ok<PFBNS>
    J = qR(edof(e,:),1:dimension)'*dNr';    %#ok<PFBNS>
    
    % Run through all Gauss points
    for k = 1:numberOfGausspoints
        indx = dimension*k-(dimension-1):dimension*k;
        detJ = det(J(:,indx)');
        if detJ < 10*eps
            error('Jacobi determinant equal or less than zero.')
        end
        dNX = (J(:,indx)')\dNr(indx,:);  %material config.
        
        %deformation gradient
        FN05 = I;
        FN05(1:dimension,1:dimension) = edN05*dNX.';
        FN1 = I;
        FN1(1:dimension,1:dimension) = edN1*dNX.';
        
        % strain energy, constitutive stresses, material tangent
        %--------------------------------------------------------------
        [~,TauN05,TauN05_v,cMatN05,errMat] = hyperelasticTB(materialObject,mapVoigtObject,FN05);
        switch errMat
            case 0 %no problem
            case 1; stretchOutOfRange=true;
        end
        WpotN1 = hyperelasticTB(materialObject,mapVoigtObject,FN1);
        strainEnergy = strainEnergy + WpotN1*detJ*gaussWeight(k);        %#ok<PFBNS>

        % residual and tangent (standard element routine)
        %--------------------------------------------------------------
        %spatial shape functions
        dNx = FN05(1:dimension,1:dimension).'\dNX;
        
        %bmatrix
        bMat = BMatrix(dNx,'mapVoigtObject',mapVoigtObject);
        
        %Residual
        Re = Re + bMat.'*TauN05_v*detJ*gaussWeight(k);
        
        %Tangent
        kgeo = dNx'*TauN05(1:dimension,1:dimension)*dNx*detJ*gaussWeight(k);
        Kgeo = zeros(numberOfDofs);
        for g = 1:dimension
            Kgeo(g:dimension:numberOfDofs,g:dimension:numberOfDofs) = kgeo;
        end
        Ke = Ke + 0.5*(bMat'*(cMatN05)*bMat)*detJ*gaussWeight(k) + Kgeo;

    end
    storeDataFE(storageFEObject,Re,Ke,globalFullEdof,e);
end
obj.ePot(setupObject.timeStep).strainEnergy = strainEnergy;
end