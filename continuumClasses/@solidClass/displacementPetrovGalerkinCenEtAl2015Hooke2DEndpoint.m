function [rData, kData, elementEnergy, array] = displacementPetrovGalerkinCenEtAl2015Hooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTPETROVGALERKINCENETAL2015HOOKE2DENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine covering linear
% mechanical processes employing a homogenous, linear-elastic, isotropic
% 'Hooke' material model (linear geometric and linear material/
% stress-strain relation).
% The formulation employs different shape functions for trial and test
% function ('PetrovGalerkin') and is suitable for static simulations.
% Another name of this element is US-ATFQ4, where US denotes the unsymmetry of
% the stiffness matrix.
%
% CALL
% displacementPetrovGalerkinCenEtAl2015Hooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which contains information like
%              time step size or plotting information.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [1, 1] for residual data.
% kData: cell-array of size [1, 1] for tangent data.
% dofs: degrees of freedom (dofs) optionally manipulated data (numerical
%       tangent)
% array: structure for storage fe-data, for more information see
%        storageFEObject.initializeArrayStress
% stressTensor: structure for storage stress tensors (postprocessing), for
%               more information see storageFEObject.initializeArrayStress
% flagNumericalTangent: flag that indicates whether the function call
%                       happens during the computation of the numerical
%                       tangent or not.
%
% REFERENCE
% https://doi.org/10.1002/nme.4899
%
% SEE ALSO
% -
%
% CREATOR(S)
% Felix Zaehringer

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;

% aquire shape functions
N_k_I = shapeFunctionObject.N_k_I;

% acquire the derivatives of the shape functions with respect to the parametric coordinates
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

% aquire the edof (element degree of freedom table)
edof = meshObject.edof;

% aquire the number of dimensions of the simulation
dimension = obj.dimension;

% number of Gausspoints
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;

% gauss points / weights
gaussPoints = shapeFunctionObject.gaussPoint;
gaussWeight = shapeFunctionObject.gaussWeight;

% aquire material data
nu = materialObject.nu;
E = materialObject.E;

% compute material matrix
switch lower(strtok(obj.materialObject.name, 'Hooke'))
    case 'esz'
        % ESZ
        C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    case 'evz'
        % EVZ
        C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
    otherwise
        error('not implemented')
end
CHat = inv(C);

x = dofs.edN1;
X = obj.qR(edof(e, :), 1:dimension)';

% compute Jacobian
JAll = computeJacobianForAllGausspoints(X, dN_xi_k_I);

% compute a-, b-, c-values
jVals = [2, 3, 4, 1];
kVals = [3, 4, 1, 2];
% a = X(1, jVals) .* X(2, kVals) - X(1, kVals) .* X(2, jVals);
b = X(2, jVals) - X(2, kVals);
c = X(1, kVals) - X(1, jVals);
% aBar = [a(3)-a(1), a(4)-a(2)];
bBar = [b(3) - b(1), b(4) - b(2)];
cBar = [c(3) - c(1), c(4) - c(2)];

% compute element area A of current element
d1 = X(:, 3) - X(:, 1);
d2 = X(:, 4) - X(:, 2);
d1Length = norm(d1);
d2Length = norm(d2);
theta = acos(d1.'*d2/(d1Length * d2Length));
A = 1 / 2 * d1Length * d2Length * sin(theta);

% compute QACM-II coordinates
nodesInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, 4);
gaussPointsInQACMIICoordinates = computeQACMIICoordinates(gaussPoints, X);
nodalPointsInQACMIICoordinates = computeQACMIICoordinates(nodesInLocalCoordinates, X);

% compute dHat matrix
dHat = computeDHatMatrix(X, nodalPointsInQACMIICoordinates, CHat, A, bBar, cBar);

% displacement vector
uN1 = x(:) - X(:);

% initialize residual
Re = rData{1};

% initialize tangent
Ke = kData{1, 1};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_k_I = computedN_X_I(dN_xi_k_I, J, k);

    % compute nodal operator matrix
    B = BMatrix(dN_X_k_I);

    % compute S and T coordinate values for current gauss point
    S = gaussPointsInQACMIICoordinates(1, k);
    T = gaussPointsInQACMIICoordinates(2, k);

    % compute 7th basic analytical solution values
    % sigmaX7 = 6/A^2*cBar(1)^2*S;
    % sigmaY7 = 6/A^2*bBar(1)^2*S;
    % sigmaXY7 = -6/A^2*bBar(1)*cBar(1)*S;
    epsX7 = 6 / (A^2) * ((cBar(1)^2 * CHat(1, 1) + bBar(1)^2 * CHat(1, 2) - bBar(1) * cBar(1) * CHat(1, 3)) * S);
    epsY7 = 6 / (A^2) * ((cBar(1)^2 * CHat(2, 1) + bBar(1)^2 * CHat(2, 2) - bBar(1) * cBar(1) * CHat(2, 3)) * S);
    epsXY7 = 6 / (A^2) * ((cBar(1)^2 * CHat(3, 1) + bBar(1)^2 * CHat(3, 2) - bBar(1) * cBar(1) * CHat(3, 3)) * S);

    % compute 8th basic analytical solution values
    % sigmaX8 = 6/A^2*cBar(2)^2*T;
    % sigmaY8 = 6/A^2*bBar(2)^2*T;
    % sigmaXY8 = -6/A^2*bBar(2)*cBar(2)*T;
    epsX8 = 6 / (A^2) * ((cBar(2)^2 * CHat(1, 1) + bBar(2)^2 * CHat(1, 2) - bBar(2) * cBar(2) * CHat(1, 3)) * T);
    epsY8 = 6 / (A^2) * ((cBar(2)^2 * CHat(2, 1) + bBar(2)^2 * CHat(2, 2) - bBar(2) * cBar(2) * CHat(2, 3)) * T);
    epsXY8 = 6 / (A^2) * ((cBar(2)^2 * CHat(3, 1) + bBar(2)^2 * CHat(3, 2) - bBar(2) * cBar(2) * CHat(3, 3)) * T);

    % PTilde matrix
    PTilde = [0, 0, 1, 0, 0, 0, epsX7, epsX8; 0, 0, 0, 0, 0, 1, epsY7, epsY8; 0, 0, 0, 1, 1, 0, epsXY7, epsXY8];

    % L matrix
    L = PTilde / dHat;

    if ~computePostData
        % RESIDUAL
        Re = Re + (B' * C * L) * uN1 * detJ * gaussWeight(k);

        % TANGENT
        % Note: Corresponds to the elements stiffness matrix
        Ke = Ke + B' * C * L * detJ * gaussWeight(k);

        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (L * uN1)' * C * (L * uN1) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        sigmaN1_v = C * L * uN1;
        stressTensor.Cauchy = voigtToMatrix(sigmaN1_v, 'stress');
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = Re;

    % pass tangent
    kData{1, 1} = Ke;
end
end

function dHat = computeDHatMatrix(nodalPointsInPhysicalCoordinates, nodalPointsInQACMIICoordinates, CHat, A, bBar, cBar)
X = nodalPointsInPhysicalCoordinates(1, :);
Y = nodalPointsInPhysicalCoordinates(2, :);

S = nodalPointsInQACMIICoordinates(1, :);
T = nodalPointsInQACMIICoordinates(2, :);

U7 = 3 / (16 * A^3) * ((cBar(1)^2 * cBar(2) * (4 * A - bBar(2) * cBar(1)) * CHat(1, 1) - bBar(1)^3 * bBar(2)^2 * CHat(2, 2) + 16 * bBar(1) * A^2 * CHat(1, 2) - bBar(1) * bBar(2)^2 * cBar(1)^2 * (CHat(1, 2) + CHat(2, 1) + CHat(3, 3)) - 16 * cBar(1) * A^2 * CHat(1, 3) + bBar(2)^2 * cBar(1)^3 * (CHat(1, 3) + CHat(3, 1)) + bBar(1)^2 * bBar(2)^2 * cBar(1) * (CHat(2, 3) + CHat(3, 2))) * S.^2 ...
    +(2 * bBar(2) * cBar(1)^4 * CHat(1, 1) + 2 * bBar(1)^4 * bBar(2) * CHat(2, 2) + 2 * bBar(1)^2 * bBar(2) * cBar(1)^2 * (CHat(1, 2) + CHat(2, 1) + CHat(3, 3)) - 2 * bBar(1) * bBar(2) * cBar(1)^3 * (CHat(1, 3) + CHat(3, 1)) - 2 * bBar(1)^3 * bBar(2) * cBar(1) * (CHat(2, 3) + CHat(3, 2))) * S .* T ...
    +(-bBar(1) * cBar(1)^4 * CHat(1, 1) - bBar(1)^5 * CHat(2, 2) - bBar(1)^3 * cBar(1)^2 * (CHat(1, 2) + CHat(2, 1) + CHat(3, 3)) + bBar(1)^2 * cBar(1)^3 * (CHat(1, 3) + CHat(3, 1)) + bBar(1)^4 * cBar(1) * (CHat(2, 3) + CHat(3, 2))) * T.^2);

V7 = 3 / (16 * A^3) * ((-cBar(1)^3 * cBar(2)^2 * CHat(1, 1) - bBar(1)^2 * bBar(2) * (4 * A + bBar(1) * cBar(2)) * CHat(2, 2) + 16 * cBar(1) * A^2 * CHat(2, 1) - bBar(1)^2 * cBar(1) * cBar(2)^2 * (CHat(1, 2) + CHat(2, 1) + CHat(3, 3)) - 16 * bBar(1) * A^2 * CHat(2, 3) + bBar(1)^3 * cBar(2)^2 * (CHat(2, 3) + CHat(3, 2)) + bBar(1) * cBar(1)^2 * cBar(2)^2 * (CHat(1, 3) + CHat(3, 1))) * S.^2 ...
    +(2 * cBar(1)^4 * cBar(2) * CHat(1, 1) + 2 * bBar(1)^4 * cBar(2) * CHat(2, 2) + 2 * bBar(1)^2 * cBar(1)^2 * cBar(2) * (CHat(1, 2) + CHat(2, 1) + CHat(3, 3)) - 2 * bBar(1) * cBar(1)^3 * cBar(2) * (CHat(1, 3) + CHat(3, 1)) - 2 * bBar(1)^3 * cBar(1) * cBar(2) * (CHat(2, 3) + CHat(3, 2))) * S .* T ...
    +(-cBar(1)^5 * CHat(1, 1) - bBar(1)^4 * cBar(1) * CHat(2, 2) - bBar(1)^2 * cBar(1)^3 * (CHat(1, 2) + CHat(2, 1) + CHat(3, 3)) + bBar(1) * cBar(1)^4 * (CHat(1, 3) + CHat(3, 1)) + bBar(1)^3 * cBar(1)^2 * (CHat(2, 3) + CHat(3, 2))) * T.^2);

U8 = 3 / (16 * A^3) * ((-cBar(2)^2 * cBar(1) * (4 * A + bBar(1) * cBar(2)) * CHat(1, 1) - bBar(2)^3 * bBar(1)^2 * CHat(2, 2) + 16 * bBar(2) * A^2 * CHat(1, 2) - bBar(2) * bBar(1)^2 * cBar(2)^2 * (CHat(1, 2) + CHat(2, 1) + CHat(3, 3)) - 16 * cBar(2) * A^2 * CHat(1, 3) + bBar(1)^2 * cBar(2)^3 * (CHat(1, 3) + CHat(3, 1)) + bBar(1)^2 * bBar(2)^2 * cBar(2) * (CHat(2, 3) + CHat(3, 2))) * T.^2 ...
    +(2 * bBar(1) * cBar(2)^4 * CHat(1, 1) + 2 * bBar(2)^4 * bBar(1) * CHat(2, 2) + 2 * bBar(2)^2 * bBar(1) * cBar(2)^2 * (CHat(1, 2) + CHat(2, 1) + CHat(3, 3)) - 2 * bBar(1) * bBar(2) * cBar(2)^3 * (CHat(1, 3) + CHat(3, 1)) - 2 * bBar(2)^3 * bBar(1) * cBar(2) * (CHat(2, 3) + CHat(3, 2))) * S .* T ...
    +(-bBar(2) * cBar(2)^4 * CHat(1, 1) - bBar(2)^5 * CHat(2, 2) - bBar(2)^3 * cBar(2)^2 * (CHat(1, 2) + CHat(2, 1) + CHat(3, 3)) + bBar(2)^2 * cBar(2)^3 * (CHat(1, 3) + CHat(3, 1)) + bBar(2)^4 * cBar(2) * (CHat(2, 3) + CHat(3, 2))) * S.^2);

V8 = 3 / (16 * A^3) * ((-cBar(1)^2 * cBar(2)^3 * CHat(1, 1) + bBar(1) * bBar(2)^2 * (4 * A - bBar(2) * cBar(1)) * CHat(2, 2) + 16 * cBar(2) * A^2 * CHat(2, 1) - bBar(2)^2 * cBar(1)^2 * cBar(2) * (CHat(1, 2) + CHat(2, 1) + CHat(3, 3)) - 16 * bBar(2) * A^2 * CHat(2, 3) + bBar(2)^3 * cBar(1)^2 * (CHat(2, 3) + CHat(3, 2)) + bBar(2) * cBar(1)^2 * cBar(2)^2 * (CHat(1, 3) + CHat(3, 1))) * T.^2 ...
    +(2 * cBar(1) * cBar(2)^4 * CHat(1, 1) + 2 * bBar(2)^4 * cBar(1) * CHat(2, 2) + 2 * bBar(2)^2 * cBar(1) * cBar(2)^2 * (CHat(1, 2) + CHat(2, 1) + CHat(3, 3)) - 2 * bBar(2) * cBar(1) * cBar(2)^3 * (CHat(1, 3) + CHat(3, 1)) - 2 * bBar(2)^3 * cBar(1) * cBar(2) * (CHat(2, 3) + CHat(3, 2))) * S .* T ...
    +(-cBar(2)^5 * CHat(1, 1) - bBar(2)^4 * cBar(2) * CHat(2, 2) - bBar(2)^2 * cBar(2)^3 * (CHat(1, 2) + CHat(2, 1) + CHat(3, 3)) + bBar(2) * cBar(2)^4 * (CHat(1, 3) + CHat(3, 1)) + bBar(2)^3 * cBar(2)^2 * (CHat(2, 3) + CHat(3, 2))) * S.^2);

dHat = zeros(8, 8);
dHat(1:2:end, :) = [ones(4, 1), zeros(4, 1), X.', zeros(4, 1), Y.', zeros(4, 1), U7.', U8.'];
dHat(2:2:end, :) = [zeros(4, 1), ones(4, 1), zeros(4, 1), X.', zeros(4, 1), Y.', V7.', V8.'];

end
