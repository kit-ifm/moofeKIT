function [rData, kData, elementEnergy, array] = easPetrovGalerkinFullIncompatibleHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
% EASPETROVGALERKINFULLINCOMPATIBLEHOOKEENDPOINT Element routine of class plateClass.


%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;
mixedFEObject = obj.mixedFEObject;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_I = shapeFunctionObject.dN0_xi_I;

numberOfNodes = size(N_k_I, 2);
gaussPoints = shapeFunctionObject.gaussPoint;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);

dimension = obj.dimension;

% aquire material data
E = materialObject.E;
nu = materialObject.nu;

% plate thickness
h = obj.h;

% nodal dofs
qN1 = dofs.edN1;
alphaN1 = dofs.edAlphaN1;

% nodal positions
ed = meshObject.nodes(edof, :).';

% compute Jacobian matrices
J0 = ed * dN0_xi_I';
detJ0 = det(J0);
J011 = J0(1, 1);
J012 = J0(1, 2);
J021 = J0(2, 1);
J022 = J0(2, 2);

% compute material matrices
% bending component
Eb = zeros(3, 3);
Eb(1, 1) = 1;
Eb(1, 2) = nu;
Eb(2, 1) = nu;
Eb(2, 2) = 1;
Eb(3, 3) = (1 - nu) / 2;
Eb = E * h^3 / (12 * (1 - nu^2)) * Eb;
% shear component
Es = eye(2, 2);
Es = 5 / 6 * E / (2 * (1 + nu)) * h * Es;

C = [Eb, zeros(3, 2); zeros(2, 3), Es];

% initialize energy
elementEnergy.strainEnergy = 0;

% compute additional shape functions
nodesInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, 4);
gaussPointsInLocalCoordinates = shapeFunctionObject.gaussPoint;
nodesInSkewCoordinates = computeSkewCoordinates(nodesInLocalCoordinates, ed, J0, shapeFunctionObject);
gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPointsInLocalCoordinates, ed, J0, shapeFunctionObject);
[M_k_I, dMr] = computeMetricShapeFunctions(obj, dimension, nodesInSkewCoordinates, gaussPointsInSkewCoordinates);
[MTilde, dMTilder] = computeShapeFunctionsTrialFunctionEnhancedStrain(nodesInSkewCoordinates, gaussPointsInSkewCoordinates, M_k_I, dMr);

% F0 matrix
F02D = F0Matrix(2, J0);
F0 = [F02D, zeros(3, 2); zeros(2, 3), J0];

% initialize residual & tangent
RX = rData{1, 1};
RA = rData{2, 1};
KXX = kData{1, 1};
KXA = kData{1, 2};
KAX = kData{2, 1};
KAA = kData{2, 2};

% compute Jacobian
JAll = computeJacobianForAllGausspoints(ed, dN_xi_k_I);

xVal = ed(1, :).';
yVal = ed(2, :).';
a1 = 1/4*[-1, 1, 1, -1].';
a2 = 1/4*[-1, -1, 1, 1].';
h = 1/4*[1, -1, 1, -1].';
C1 = 1/detJ0*(a2.'*yVal*h.'*xVal-a2.'*xVal*h.'*yVal);
C2 = 1/detJ0*(-a1.'*yVal*h.'*xVal+a1.'*xVal*h.'*yVal);

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % derivatives with respect to the physical coordinates
    indx = dimension * k - (dimension - 1):dimension * k;
    dMx = J0' \ dMr(indx, :);
    dMTildex = J0' \ dMTilder(indx, :);

    xi = gaussPoints(1, k);
    eta = gaussPoints(2, k);
    Er = zeros(5, 6);
    Er(1, 3) = xi;
    Er(2, 6) = eta;
    Er(3, 4) = eta;
    Er(3, 5) = xi;
    Er(4, 1) = xi;
    Er(5, 2) = eta;
    Er(4, 5) = xi*eta-1/3*C2*eta;
    Er(5, 4) = xi*eta-1/3*C1*xi;

    % shape functions for the enhanced part of the strain field
    Hr = zeros(5, 6);
    % bending part
    
    % shear part
    Hr(4:5, 1:2) = dMTilder(indx, 1:2);
    Hr(4, 3) = -(J012*dMTilder(indx(1), 3)+J011*dMTilder(indx(1), 5))+J011*MTilde(k, 1);
    Hr(5, 3) = -(J012*dMTilder(indx(2), 3)+J011*dMTilder(indx(2), 5))+J012*MTilde(k, 1);
    Hr(4, 4) = -(J011*dMTilder(indx(1), 4)+J012*dMTilder(indx(1), 6))+J011*MTilde(k, 2);
    Hr(5, 4) = -(J011*dMTilder(indx(2), 4)+J012*dMTilder(indx(2), 6))+J012*MTilde(k, 2);
    Hr(4, 5) = -(J022*dMTilder(indx(1), 3)+J021*dMTilder(indx(1), 5))+J021*MTilde(k, 1);
    Hr(5, 5) = -(J022*dMTilder(indx(2), 3)+J021*dMTilder(indx(2), 5))+J022*MTilde(k, 1);
    Hr(4, 6) = -(J021*dMTilder(indx(1), 4)+J022*dMTilder(indx(1), 6))+J021*MTilde(k, 2);
    Hr(5, 6) = -(J021*dMTilder(indx(2), 4)+J022*dMTilder(indx(2), 6))+J022*MTilde(k, 2);

    % B-Matrix
    % bending component
    Bb = zeros(3, 3*numberOfNodes);
    Bb(1, 2:3:end) = dN_X_I(1, :);
    Bb(2, 3:3:end) = dN_X_I(2, :);
    Bb(3, 2:3:end) = dN_X_I(2, :);
    Bb(3, 3:3:end) = dN_X_I(1, :);

    Lb = zeros(3, 3*numberOfNodes);
    Lb(1, 2:3:end) = dMx(1, :);
    Lb(2, 3:3:end) = dMx(2, :);
    Lb(3, 2:3:end) = dMx(2, :);
    Lb(3, 3:3:end) = dMx(1, :);
    
    % shear component
    Bs=zeros(2,3*numberOfNodes);
    Bs(1, 1:3:end) = dN_X_I(1, :);
    Bs(2, 1:3:end) = dN_X_I(2, :);
    Bs(1, 2:3:end) = N_k_I(k, :);
    Bs(2, 3:3:end) = N_k_I(k, :);

    Ls=zeros(2,3*numberOfNodes);
    Ls(1, 1:3:end) = dMx(1, :);
    Ls(2, 1:3:end) = dMx(2, :);
    Ls(1, 2:3:end) = M_k_I(k, :);
    Ls(2, 3:3:end) = M_k_I(k, :);

    B = [Bb; Bs];
    L = [Lb; Ls];

    G = detJ0 / detJ * (F0.' \ Er);


    Hold = F0.' \ Hr;

    H = zeros(5, 6);
    H(1, 3) = dMTildex(1, 1);
    H(1, 4) = dMTildex(1, 2);
    H(2, 5) = dMTildex(2, 1);
    H(2, 6) = dMTildex(2, 2);
    H(3, 3) = dMTildex(2, 1);
    H(3, 4) = dMTildex(2, 2);
    H(3, 5) = dMTildex(1, 1);
    H(3, 6) = dMTildex(1, 2);
    H(4, 1) = dMTildex(1, 1);
    H(4, 2) = dMTildex(1, 2);
    H(4, 3) = -(J012*dMTildex(1, 3)+J011*dMTildex(1, 5))+MTilde(k, 1);
    H(4, 4) = -(J011*dMTildex(1, 4)+J012*dMTildex(1, 6))+MTilde(k, 2);
    H(4, 5) = -(J022*dMTildex(1, 3)+J021*dMTildex(1, 5));
    H(4, 6) = -(J021*dMTildex(1, 4)+J022*dMTildex(1, 6));
    H(5, 1) = dMTildex(2, 1);
    H(5, 2) = dMTildex(2, 2);
    H(5, 3) = -(J012*dMTildex(2, 3)+J011*dMTildex(2, 5));
    H(5, 4) = -(J011*dMTildex(2, 4)+J012*dMTildex(2, 6));
    H(5, 5) = -(J022*dMTildex(2, 3)+J021*dMTildex(2, 5))+MTilde(k, 1);
    H(5, 6) = -(J021*dMTildex(2, 4)+J022*dMTildex(2, 6))+MTilde(k, 2);
    
    if ~computePostData
        % Tangent
        KXX = KXX + B.' * C * L * detJ * gaussWeight(k);
        KXA = KXA + B.' * C * H * detJ * gaussWeight(k);
        KAX = KAX + G.' * C * L * detJ * gaussWeight(k);
        KAA = KAA + G.' * C * H * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        eps = L * qN1(:) + H * alphaN1(:);

        m = C * eps;
        stressTensor.Cauchy = [m(1), m(3), m(5); m(3), m(2), m(4); m(5), m(4), 0];
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension+1);
    end
end

%% RESIDUAL
RX = RX + KXX * qN1(:) + KXA * alphaN1(:);
RA = RA + KAX * qN1(:) + KAA * alphaN1(:);

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = RX;
    rData{2} = RA;
    kData{1, 1} = KXX;
    kData{1, 2} = KXA;
    kData{2, 1} = KAX;
    kData{2, 2} = KAA;
end
end

function [MTilde, dMTilder] = computeShapeFunctionsTrialFunctionEnhancedStrain(nodesInSkewCoordinates, gaussPointsInSkewCoordinates, M, dMr)
% Computes the shape functions for the trial function for the enhanced part
% of the strain field
MTilde = zeros(size(gaussPointsInSkewCoordinates, 2), 6);
MTilde(:, 1:2) = (gaussPointsInSkewCoordinates.^2).' - M * (nodesInSkewCoordinates.^2).';
MTilde(:, 3) = (gaussPointsInSkewCoordinates(1, :).^2.*gaussPointsInSkewCoordinates(2, :)).' - M * (nodesInSkewCoordinates(1, :).^2.*nodesInSkewCoordinates(2, :)).';
MTilde(:, 4) = (gaussPointsInSkewCoordinates(1, :).*gaussPointsInSkewCoordinates(2, :).^2).' - M * (nodesInSkewCoordinates(1, :).*nodesInSkewCoordinates(2, :).^2).';
MTilde(:, 5) = 1/3 * (gaussPointsInSkewCoordinates(1, :).^3).' - 1/3 * M * (nodesInSkewCoordinates(1, :).^3).';
MTilde(:, 6) = 1/3 * (gaussPointsInSkewCoordinates(2, :).^3).' - 1/3 * M * (nodesInSkewCoordinates(2, :).^3).';

dMTilder = zeros(2*size(gaussPointsInSkewCoordinates, 2), 6);
dMTilder(1:2:end, 1) = 2*gaussPointsInSkewCoordinates(1, :).' - dMr(1:2:end, :) * (nodesInSkewCoordinates(1, :).^2).';
dMTilder(2:2:end, 1) = - dMr(2:2:end, :) * (nodesInSkewCoordinates(1, :).^2).';
dMTilder(1:2:end, 2) = - dMr(1:2:end, :) * (nodesInSkewCoordinates(2, :).^2).';
dMTilder(2:2:end, 2) = 2*gaussPointsInSkewCoordinates(2, :).' - dMr(2:2:end, :) * (nodesInSkewCoordinates(2, :).^2).';
dMTilder(1:2:end, 3) = 2*(gaussPointsInSkewCoordinates(1, :).*gaussPointsInSkewCoordinates(2, :)).' - dMr(1:2:end, :) * (nodesInSkewCoordinates(1, :).^2.*nodesInSkewCoordinates(2, :)).';
dMTilder(2:2:end, 3) = (gaussPointsInSkewCoordinates(1, :).^2).' - dMr(2:2:end, :) * (nodesInSkewCoordinates(1, :).^2.*nodesInSkewCoordinates(2, :)).';
dMTilder(1:2:end, 4) = (gaussPointsInSkewCoordinates(2, :).^2).' - dMr(1:2:end, :) * (nodesInSkewCoordinates(1, :).*nodesInSkewCoordinates(2, :).^2).';
dMTilder(2:2:end, 4) = 2*(gaussPointsInSkewCoordinates(1, :).*gaussPointsInSkewCoordinates(2, :)).' - dMr(2:2:end, :) * (nodesInSkewCoordinates(1, :).*nodesInSkewCoordinates(2, :).^2).';
dMTilder(1:2:end, 5) = (gaussPointsInSkewCoordinates(1, :).^2).' - 1/3 * dMr(1:2:end, :) * (nodesInSkewCoordinates(1, :).^3).';
dMTilder(2:2:end, 6) =(gaussPointsInSkewCoordinates(2, :).^2).' - 1/3 * dMr(2:2:end, :) * (nodesInSkewCoordinates(2, :).^3).';
end
