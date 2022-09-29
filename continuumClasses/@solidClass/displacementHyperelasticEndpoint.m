function [rData, kData, elementEnergy, array] = displacementHyperelasticEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% 21.02.2022:
% This formulation does not work currently (missing function
% 'hyperelasticTB')!

% out = disp_logStrain_endPoint(obj,'PropertyName',PropertyValue)
%
% Description:
% -various hyperelastic laws,
% -evaluated at time n+1, i.e. implicid euler method.
% -spatial formulation in Tau and left Cauchy Green
%
% 19.4.2021 Robin Pfefferkorn

%% SETUP
% load objects
meshObject = obj.meshObject;
shapeFunctionObject = obj.shapeFunctionObject;
storageFEObject = obj.storageFEObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
mixedFEObject = obj.mixedFEObject;

% aquire general data
N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;

edof = meshObject.edof;
numberOfDofs = size(meshObject.globalFullEdof, 2);

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

dimension = obj.dimension;

selectMapVoigt(mapVoigtObject, dimension, 'symmetric');

I = eye(3);

% aquire the nodal values of the variables for the current element
edRef = obj.qR(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;

% compute Jacobian matrix
J = edRef * dNr';

% initialize
Re = rData{1};
Ke = kData{1};
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    indx = dimension * k - (dimension - 1):dimension * k;
    detJ = det(J(:, indx)');
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNX = (J(:, indx)') \ dNr(indx, :);

    % deformation gradient
    FN1 = I;
    FN1(1:dimension, 1:dimension) = edN1 * dNX.';

    % strain energy, constitutive stresses, material tangent
    [Wpot, TauN1, TauN1_v, cMat, errMat] = hyperelasticTB(materialObject, mapVoigtObject, FN1);
    
    if ~computePostData
        % spatial shape functions
        dNx = FN1(1:dimension, 1:dimension).' \ dNX;

        % B-matrix
        bMat = BMatrix(dNx, 'mapVoigtObject', mapVoigtObject);

        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + Wpot * detJ * gaussWeight(k);

        % RESIDUAL
        Re = Re + bMat.' * TauN1_v * detJ * gaussWeight(k);

        % TANGENT
        kgeo = dNx' * TauN1(1:dimension, 1:dimension) * dNx * detJ * gaussWeight(k);
        Kgeo = zeros(numberOfDofs);
        for g = 1:dimension
            Kgeo(g:dimension:numberOfDofs, g:dimension:numberOfDofs) = kgeo;
        end
        Ke = Ke + (bMat' * (cMat) * bMat) * detJ * gaussWeight(k) + Kgeo;
    else
        % STRESS COMPUTATION
        sigma = det(FN1)^-1 * TauN1;
        tempSpann = selectStress(sigma, computeStresses, 3);
        se = se + N(k, :)' * tempSpann * detJ * gaussWeight(k);
        Me = Me + (N(k, :)' * N(k, :)) * detJ * gaussWeight(k);
    end
end

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = Re;
    kData{1} = Ke;
end
end