function [rData, kData, elementEnergy, array] = easKasperTaylorHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% EASKASPERTAYLORHOOKEENDPOINT Element routine of class axisymmetricSolidClass.
%
% FORMULATION
% This is a 'eas'-based finite element routine covering linear
% mechanical processes employing a elastic, isotropic Hooke
% ('Hooke') model.
% The routine is suitable for static and dynamic simulations where for the
% latter the backward Euler integration scheme is used ('Endpoint').
%
% CALL
% easKasperTaylorHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type axisymmetricSolidClass,
%      e.g. axisymmetricSolidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [2, 1] for residual data of
%        every field, here: (X, alpha)
% kData: cell-array of size [2, 2] for
%        tangent data of every field, here: (X, alpha)
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
% https://doi.org/10.1002/nme.373
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

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_I = shapeFunctionObject.dN0_xi_I;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;
gaussPoints = shapeFunctionObject.gaussPoint;

edof = meshObject.edof;

dimension = obj.dimension;

% aquire material data
E = materialObject.E;
nu = materialObject.nu;
C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu,  nu, 0; nu, 1 - nu, nu, 0; nu, nu, 1 - nu, 0; 0, 0, 0, (1 - 2 * nu) / 2];

% aquire the nodal values of the variables for the current element
edR = obj.qR(edof(e, :), 1:dimension).';
edN1 = dofs.edN1;
alphaN1e = dofs.edAlphaN1.';
uN1 = edN1(:) - edR(:);

% compute Jacobian matrices
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);
J0 = edR * dN0_xi_I.';

% initialize tangent
KDD = kData{1, 1};
KDA = kData{1, 2};
KAD = kData{2, 1};
KAA = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

% F0 Matrix
F0([1, 2, 4], [1, 2, 4]) = F0Matrix(dimension, J0);
F0(3, 3) = 1;

% element specific stuff
T = eye(3, 3);
T(1:2, 1:2) = J0;

r1 = edR(1, 1);
r2 = edR(1, 2);
r3 = edR(1, 3);
r4 = edR(1, 4);
rBar = r1+r2+r3+r4;

xiBar = (r2-r1+r3-r4)/(3*rBar);
etaBar = (r4-r1+r3-r2)/(3*rBar);
lambda1 = 2/(9*rBar)*((r1+r2)^2+(r3+r4)^2+4*(r1+r2)*(r3+r4));
lambda2 = 2/(9*rBar)*((r1+r4)^2+(r2+r3)^2+4*(r1+r4)*(r2+r3));
lambda3 = 4/(9*rBar)*(r1*r3-r2*r4);


%% GAUSS LOOP
BBar = zeros(4, 8);
V = 0;
JAvg = zeros(dimension);
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    r = N_k_I(k, :) * edR(1, :).';

    B = zeros(4, 8);
    B([1,2,4], :) = BMatrix(dN_X_I);
    B(3, 1:2:end) = 1 / r * N_k_I(k, :);

    BBar = BBar + 2 * pi * B * r * detJ * gaussWeight(k);
    V = V + 2 * pi * r * detJ * gaussWeight(k);
    JAvg = JAvg + 2 * pi * J * r * detJ * gaussWeight(k);
end
BBar = 1/V*BBar;
% 
% JAvg = 1/V*JAvg;
% T = eye(3, 3);
% T(1:2, 1:2) = JAvg';
% t = inv(T);

H = zeros(4, 8);
HAlt = zeros(4, 8);
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    xi = gaussPoints(1, k);
    eta = gaussPoints(2, k);

    r = N_k_I(k, :) * edR(1, :).';

    Gamma1 = 1/(lambda3^2-lambda1*lambda2)*((xi-xiBar)*lambda1-(eta-etaBar)*lambda3);
    Gamma2 = 1/(lambda3^2-lambda1*lambda2)*((eta-etaBar)*lambda2-(xi-xiBar)*lambda3);

    B = zeros(4, 8);
    B([1,2,4], :) = BMatrix(dN_X_I);
    B(3, 1:2:end) = 1 / r * N_k_I(k, :);

    BStar = B-BBar;

    Htemp = zeros(4, 8);
    Htemp(1,1:2:end) = (eta-etaBar)/lambda1*(T(1,1)^2*BStar(1, 1:2:end)+T(1, 1)*T(2, 1)*BStar(2, 2:2:end));
    Htemp(1,2:2:end) = (eta-etaBar)/lambda1*(T(1,1)*T(2, 1)*BStar(1, 1:2:end)+T(2, 1)^2*BStar(2, 2:2:end));
    Htemp(2,1:2:end) = (xi-xiBar)/lambda2*(T(1,2)^2*BStar(1, 1:2:end)+T(1, 2)*T(2, 2)*BStar(2, 2:2:end));
    Htemp(2,2:2:end) = (xi-xiBar)/lambda2*(T(1,2)*T(2, 2)*BStar(1, 1:2:end)+T(2, 2)^2*BStar(2, 2:2:end));
%     Htemp(3,1:2:end) = - Gamma1*N_k_I(k,:)/r;
%     Htemp(4,1:2:end) = - Gamma2*N_k_I(k,:)/r;
    Htemp(3,1:2:end) = - Gamma1*BStar(3, 1:2:end);
    Htemp(4,1:2:end) = - Gamma2*BStar(3, 1:2:end);

    H = H + 2 * pi * Htemp * r * detJ * gaussWeight(k);
    
    E1 = zeros(4, 4);
    E1(1, 1) = eta-etaBar;
    E1(2, 2) = xi-xiBar;
    E1(3, 3) = xi-xiBar;
    E1(3, 4) = eta-etaBar;
    
    HGamma = [lambda1, 0, 0, 0; 0, lambda2, 0, 0; 0, 0, lambda2, lambda3; 0, 0, lambda3, lambda1];
    HAlt = HAlt + 2*pi* (inv(HGamma) * E1.' * F0.' * BStar) * r * detJ * gaussWeight(k);
%     assert(max(max(H-HAlt))<1e-15);
end
for k = 1:numberOfGausspoints
    [~, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);

    xi = gaussPoints(1, k);
    eta = gaussPoints(2, k);

    r = N_k_I(k, :) * edR(1, :).';
    
    Bu3 = zeros(4, 4);
    Bu3(1, 1) = eta-etaBar;
    Bu3(2, 2) = xi-xiBar;
    Bu3(3, 3) = xi-xiBar;
    Bu3(3, 4) = eta-etaBar;

    Bu4 = H/(2*pi);
    
    Bu5 = [(xi-xiBar)-lambda3/lambda1*(eta-etaBar), 0; 0, (eta-etaBar)-lambda3/lambda2*(xi-xiBar); 0, 0; 0,0];

    Bu = BBar+1/detJ*(F0.'\Bu3*Bu4);
    BAlpha = 1/detJ*(F0.'\Bu5);

    if ~computePostData
        % TANGENT
        KDD = KDD + 2 * pi * (Bu.' * C * Bu) * r * detJ * gaussWeight(k);
        KDA = KDA + 2 * pi * (Bu.' * C * BAlpha) * r * detJ * gaussWeight(k);
        KAD = KAD + 2 * pi * (BAlpha.' * C * Bu) * r * detJ * gaussWeight(k);
        KAA = KAA + 2 * pi * (BAlpha.' * C * BAlpha) * r * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        sigmaN1_v = C * (Bu * uN1 + BAlpha * alphaN1e);
        stressTensor.Cauchy = [sigmaN1_v(1), sigmaN1_v(3), 0; sigmaN1_v(3), sigmaN1_v(2), 0; 0, 0, sigmaN1_v(4)];
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension+1);
    end
end

%% RESIDUAL
RD = KDD * uN1 + KDA * alphaN1e;
RA = KAD * uN1 + KAA * alphaN1e;

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = RD;
    rData{2} = RA;
    kData{1, 1} = KDD;
    kData{1, 2} = KDA;
    kData{2, 1} = KAD;
    kData{2, 2} = KAA;
end
end