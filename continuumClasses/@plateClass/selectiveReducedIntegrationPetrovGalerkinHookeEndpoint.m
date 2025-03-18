function [rData, kData, elementEnergy, array] = selectiveReducedIntegrationPetrovGalerkinHookeEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, ~)
% homogenous linear-elastic isotropic strain-energy function, evaluated at time n+1, i.e. implicid euler method.
% 08.02.2015 Marlon Franke: first esra version
% 19.08.2021 Marlon Franke: moofeKIT version
%creates the residual and the tangent of the given obj.
%

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;
mixedFEObject = obj.mixedFEObject;

% aquire general data
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_I = shapeFunctionObject.dN0_xi_I;

numberOfNodes = size(dN_xi_k_I, 3);

edof = meshObject.edof(e, :);

dimension = obj.dimension;

% aquire material data
E = materialObject.E;
nu = materialObject.nu;

% plate thickness
h = obj.h;

% nodal dofs
qN1 = dofs.edN1;

% nodal positions
ed = meshObject.nodes(edof, :).';

% Jacobian
J0 = ed * dN0_xi_I';

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

% initialize energy
elementEnergy.strainEnergy = 0;

% initialize residual & tangent
RX = rData{1, 1};
KXX = kData{1, 1};

for integrationType = 1:2
    %shape functions
    if integrationType == 1
        numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
        gaussWeight = shapeFunctionObject.gaussWeight;
        N_k_I = shapeFunctionObject.N_k_I;
        dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

        nodesInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, 4);
        gaussPointsInLocalCoordinates = shapeFunctionObject.gaussPoint;
        nodesInSkewCoordinates = computeSkewCoordinates(nodesInLocalCoordinates, ed, J0, shapeFunctionObject);
        gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPointsInLocalCoordinates, ed, J0, shapeFunctionObject);
        [M_k_I, dMr] = computeMetricShapeFunctions(obj, dimension, nodesInSkewCoordinates, gaussPointsInSkewCoordinates);
    else
        numberOfGausspoints = size(mixedFEObject.shapeFunctionObject.gaussWeight, 2);
        gaussWeight = mixedFEObject.shapeFunctionObject.gaussWeight;
        N_k_I = mixedFEObject.shapeFunctionObject.N_k_I;
        dN_xi_k_I = mixedFEObject.shapeFunctionObject.dN_xi_k_I;

%         numberOfGausspoints = 4;
%         [gaussPoints, gaussWeight] = gaussPointsAndWeights(dimension,numberOfGausspoints,obj.elementGeometryType);
%         [N_k_I, dN_xi_k_I] = computeLagrangeShapeFunction(dimension, numberOfNodes, numberOfGausspoints,gaussPoints);

        nodesInLocalCoordinates = elementNodesInLocalCoordinates(dimension, obj.elementGeometryType, 4);
        gaussPointsInLocalCoordinates = mixedFEObject.shapeFunctionObject.gaussPoint;
%         gaussPointsInLocalCoordinates = gaussPoints;
        nodesInSkewCoordinates = computeSkewCoordinates(nodesInLocalCoordinates, ed, J0, mixedFEObject.shapeFunctionObject);
        gaussPointsInSkewCoordinates = computeSkewCoordinates(gaussPointsInLocalCoordinates, ed, J0, mixedFEObject.shapeFunctionObject);
        [M_k_I, dMr] = computeMetricShapeFunctions(obj, dimension, nodesInSkewCoordinates, gaussPointsInLocalCoordinates);
    end

    %% GAUSS LOOP
    for k = 1:numberOfGausspoints
        [detJ, detJStruct, dN_X_I, ~] = computeAllJacobian(ed, ed, ed, dN_xi_k_I, k, setupObject);

        % derivatives with respect to the physical coordinates
        index = dimension * (k - 1) + 1:dimension * k;
        dM_X_k_I = J0' \ dMr(index, :);

        % B-Matrix
        % bending component
        Bb = zeros(3, 3*numberOfNodes);
        Bb(1, 2:3:end) = dN_X_I(1, :);
        Bb(2, 3:3:end) = dN_X_I(2, :);
        Bb(3, 2:3:end) = dN_X_I(2, :);
        Bb(3, 3:3:end) = dN_X_I(1, :);

        Lb = zeros(3, 3*numberOfNodes);
        Lb(1, 2:3:end) = dM_X_k_I(1, :);
        Lb(2, 3:3:end) = dM_X_k_I(2, :);
        Lb(3, 2:3:end) = dM_X_k_I(2, :);
        Lb(3, 3:3:end) = dM_X_k_I(1, :);

        % shear component
        Bs = zeros(2, 3*numberOfNodes);
        Bs(1, 1:3:end) = dN_X_I(1, :);
        Bs(2, 1:3:end) = dN_X_I(2, :);
        Bs(1, 2:3:end) = N_k_I(k, :);
        Bs(2, 3:3:end) = N_k_I(k, :);

        Ls = zeros(2, 3*numberOfNodes);
        Ls(1, 1:3:end) = dM_X_k_I(1, :);
        Ls(2, 1:3:end) = dM_X_k_I(2, :);
        Ls(1, 2:3:end) = M_k_I(k, :);
        Ls(2, 3:3:end) = M_k_I(k, :);

        if ~computePostData
            % Tangent
            % bending component
            if integrationType == 1
                Kb = Bb' * Eb * Lb * detJ * gaussWeight(k);
                KXX = KXX + Kb;
            elseif integrationType == 2
                % shear component
                Ks = Bs' * Es * Ls * detJ * gaussWeight(k);
                KXX = KXX + Ks;
            end

        else
            % stress at gausspoint
            kappa = Lb * qN1(:);
            gamma = Ls * qN1(:);
    
            m = Eb * kappa;
            q = Es * gamma;
            stressTensor.Cauchy = [m(1), m(3), q(1); m(3), m(2), q(2); q(1), q(2), 0];
            array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension+1);
        end
    end

end

%% RESIDUAL
RX = KXX * qN1(:);

%% PASS COMPUTATION DATA
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end