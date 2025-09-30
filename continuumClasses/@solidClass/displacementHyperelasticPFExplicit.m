function displacementHyperelasticPFExplicit(obj,setupObject,varargin)
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
numberOfDofs = size(globalFullEdof,2);
materialObject = obj.materialObject;
selectMapVoigt(obj,'unsymmetric');

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

    %dofs and jacobian
    edN = qN(edof(e,:),1:dimension).'; %#ok<PFBNS>
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
        FN = I;
        FN(1:dimension,1:dimension) = edN*dNX.';

        % strain energy, constitutive stresses, material tangent
        %--------------------------------------------------------------
        [WpotN,~,PN_v,~,errMat] = hyperelasticPF(materialObject,obj,FN);
        switch errMat
            case 0 %no problem
            case 1; stretchOutOfRange=true;
        end
        strainEnergy = strainEnergy + WpotN*detJ*gaussWeight(k);        %#ok<PFBNS>

        % residual and tangent (standard element routine)
        %--------------------------------------------------------------

        %bmatrix
        BMat = BMatrix(dNX,'type',obj.mapType);

        %Residual
        Re = Re + BMat.'*PN_v*detJ*gaussWeight(k);

    end
    storeDataFE(obj,Re,Ke,globalFullEdof,e);
end
obj.elementData(setupObject.timeStep).strainEnergy = strainEnergy;
end