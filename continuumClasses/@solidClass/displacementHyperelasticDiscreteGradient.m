function displacementHyperelasticDiscreteGradient(obj,setupObject,varargin)
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
discreteGradientFlag = strcmpi(setupObject.integrator,'discreteGradient');

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
    edN = qN(edof(e,:),1:dimension).'; %#ok<PFBNS>
    edN1 = qN1(edof(e,:),1:dimension).'; %#ok<PFBNS>
    edN05 = 0.5*(edN + edN1); %#ok<PFBNS>
    J = qR(edof(e,:),1:dimension)'*dNr';    %#ok<PFBNS>
    
    % test for activation of discrete gradient
    activeDisGra = false;
    if discreteGradientFlag        
        for k = 1:numberOfGausspoints
            indx = dimension*k-(dimension-1):dimension*k;
            detJ = det(J(:,indx)');
            if detJ < 10*eps
                error('Jacobi determinant equal or less than zero.')
            end
            dNX = (J(:,indx)')\dNr(indx,:);
            FN = edN*dNX';
            FN1 = edN1*dNX';
            % Cauchy-Green-Tensor
            CN  = FN'*FN;
            CN1 = FN1'*FN1;
            % Frobenius norm
            if norm(CN1-CN,'fro')>1e-11
                activeDisGra = true;
                break;
            end            
        end
    end    
    
    % Run through all Gauss points
    for k = 1:numberOfGausspoints
        indx = dimension*k-(dimension-1):dimension*k;
        detJ = det(J(:,indx)');
        if detJ < 10*eps
            error('Jacobi determinant equal or less than zero.')
        end
        dNX = (J(:,indx)')\dNr(indx,:);  %material config.
        
        %deformation gradient
        FN = I;
        FN(1:dimension,1:dimension) = edN*dNX.';
        FN05 = I;
        FN05(1:dimension,1:dimension) = edN05*dNX.';
        FN1 = I;
        FN1(1:dimension,1:dimension) = edN1*dNX.';
        
        % strain energy, constitutive stresses, material tangent
        %--------------------------------------------------------------
        [WpotN1,Salg,Salg_v,CCalg,errMat] = hyperelasticDiscreteGradient(materialObject,mapVoigtObject,FN.'*FN,FN1.'*FN1,activeDisGra);
        switch errMat
            case 0 %no problem
            case 1; stretchOutOfRange=true;
        end
        strainEnergy = strainEnergy + WpotN1*detJ*gaussWeight(k);        %#ok<PFBNS>

        % residual and tangent (standard element routine)
        %--------------------------------------------------------------
        %bmatrix
        BMatN05 = BMatrix(dNX,FN05);
        BMatN1 = BMatrix(dNX,FN1);
        
        %Residual
        Re = Re + BMatN05.'*Salg_v*detJ*gaussWeight(k);
        
        %Tangent
        kgeo = dNX'*Salg(1:dimension,1:dimension)*dNX*detJ*gaussWeight(k);
        Kgeo = zeros(numberOfDofs);
        for g = 1:dimension
            Kgeo(g:dimension:numberOfDofs,g:dimension:numberOfDofs) = kgeo;
        end
        Ke = Ke + 0.5*(BMatN05'*(CCalg)*BMatN1)*detJ*gaussWeight(k) + Kgeo;

    end
    storeDataFE(storageFEObject,Re,Ke,globalFullEdof,e);
end
obj.ePot(setupObject.timeStep).strainEnergy = strainEnergy;
end