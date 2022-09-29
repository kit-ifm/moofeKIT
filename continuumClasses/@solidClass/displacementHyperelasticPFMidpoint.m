function displacementHyperelasticPFMidpoint(obj,setupObject,varargin)
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

%% get & prepare Data
%Shape functions
globalFullEdof = obj.globalFullEdof;
edof = obj.edof;
numberOfGausspoints = obj.numberOfGausspoints;
gaussWeight = obj.shapeFunctions.gaussWeight;
N = obj.shapeFunctions.N;
dNr = obj.shapeFunctions.dNr;
dimension = obj.dimension;
numberOfElements = size(globalFullEdof,1);
qR = obj.qR;
qN = obj.qN;
qN1 = obj.qN1;
qN05 = 0.5*(qN + qN1);
numberOfDofs = size(globalFullEdof,2);
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
selectMapVoigt(mapVoigtObject,dimension,'unsymmetric');

%identity tensors
I = eye(3);

%% output
strainEnergy = 0;
initializeDataFE(obj);
%% Element Loop
%---------------------------------------------
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
        [~,PN05,PN05_v,CMatN05,errMat] = hyperelasticPF(materialObject,mapVoigtObject,FN05);
        switch errMat
            case 0 %no problem
            case 1; stretchOutOfRange=true;
        end
        WpotN1 = hyperelasticPF(materialObject,mapVoigtObject,FN1);
        strainEnergy = strainEnergy + WpotN1*detJ*gaussWeight(k);        %#ok<PFBNS>

        % residual and tangent (standard element routine)
        %--------------------------------------------------------------
        %bmatrix
        BMat = BMatrix(dNX,'type',mapVoigtObject.mapType);
        
        %Residual
        Re = Re + BMat.'*PN05_v*detJ*gaussWeight(k);
        
        %Tangent
        Ke = Ke + 0.5*(BMat'*(CMatN05)*BMat)*detJ*gaussWeight(k);

    end
    storeDataFE(obj,Re,Ke,globalFullEdof,e);
end
obj.ePot(setupObject.timeStep).strainEnergy = strainEnergy;
end