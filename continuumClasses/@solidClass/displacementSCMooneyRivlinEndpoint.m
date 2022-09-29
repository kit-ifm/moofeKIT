function [rData, kData, elementEnergy, array] = displacementSCMooneyRivlinEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTSCMOONEYRIVLINENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine  covering nonlinear
% mechanical processes employing a hyperelastic, isotropic Mooney-Rivlin
% ('MooneyRivlin') model (attention: nonlinear geometric and material/
% stress-strain relation).
% The formulation is based on a polyconvexity inspired strain-energy
% density function.
% Implementation is due to work-conjugated 2nd PK-stress tensor and Cauchy-
% Green strain tensor ('SC').
% The routine is suitable for static and dynamic simulation where for the
% latter the backward Euler integration scheme is used ('Endpoint').
%
% CALL
% displacementSCMooneyRivlinEndpoint(obj,setupObject,computePostData)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which cotains informations like
%              time step size or plotting informations.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
%
% REFERENCE
% https://doi.org/10.1016/j.cma.2018.01.013
%
% SEE ALSO
% displacementSCMooneyRivlinMidpoint,
% displacementSCMooneyRivlinDiscreteGradient
%
% CREATOR(S)
% Marlon Franke

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
storageFEObject = obj.storageFEObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
%   mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
% element degree of freedom tables and more
edof = meshObject.edof;
globalFullEdof = meshObject.globalFullEdof;
numberOfElements = size(globalFullEdof, 1);
numberOfDOFs = size(globalFullEdof, 2);
dimension = obj.dimension;
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
% nodal dofs
qR = obj.qR;
qN1 = obj.qN1;
% material data and voigt notation
selectMapVoigt(mapVoigtObject, dimension, 'symmetric');
a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;
I = eye(dimension);

%% Create residual and tangent
elementEnergy.strainEnergy = 0;
RX = rData{1};
KXX = kData{1, 1};

edN1 = dofs.edN1;
J = qR(edof(e, :), 1:dimension)' * dNr';
JN1 = edN1 * dNr';
% Run through all Gauss points
for k = 1:numberOfGausspoints
    index = dimension * k - (dimension - 1):dimension * k;
    detJ = det(J(:, index)');
    detJN1 = det(JN1(:, index)');
    if detJ < 10 * eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:, index)') \ dNr(index, :);
    % Deformation gradient
    FN1 = edN1 * dNx';
    % B-matrix (current configuration)
    BN1 = BMatrix(dNx, FN1);
    % Right Cauchy-Green tensor
    CN1 = FN1' * FN1;
    % Cofactor
    GN1 = 0.5 * wedge(CN1, CN1);
    % Third invariant
    I3N1 = det(CN1);
    % Derivative of the strain energy function
    Sigma_C = a * I;
    Sigma_G = b * I;
    Sigma_I = -d / (2 * I3N1) + c / 2 * (1 - 1 / sqrt(I3N1));
    Sigma_I_I = d / (2 * I3N1^2) + c / (4 * (I3N1)^(3 / 2));
    % Second Piola Kirchhoff stress tensor
    SN1 = 2 * (Sigma_C + wedge(Sigma_G, CN1) + Sigma_I * GN1);
    SN1_v = [SN1(1, 1); SN1(2, 2); SN1(3, 3); SN1(1, 2); SN1(2, 3); SN1(1, 3)];
    if ~computePostData
        % Residual
        RX = RX + BN1' * SN1_v * detJ * gaussWeight(k);
        % Tangent
        % Derivative of wedge(Sigma_G,CN05)
        Kmat1 = secDiffOperator(Sigma_G);
        % Derivative of Sigma_I*1/3*(wedge(CN05,CN05)+GN05)  part I
        Kmat2 = Sigma_I * secDiffOperator(CN1);
        % Derivative of Sigma_I*1/3*(wedge(CN05,CN05)+GN05)  part II
        GN1_v = [GN1(1, 1); GN1(2, 2); GN1(3, 3); GN1(1, 2); GN1(2, 3); GN1(1, 3)];
        Kmat3 = Sigma_I_I * (GN1_v * GN1_v');
        % Assembly of elasticity tensor
        ELA = 4 * (Kmat1 + Kmat2 + Kmat3);
        A1 = dNx' * SN1 * dNx * detJ * gaussWeight(k);
        MAT = zeros(numberOfDOFs);
        for g = 1:dimension
            MAT(g:dimension:numberOfDOFs, g:dimension:numberOfDOFs) = A1;
        end
        KXX = KXX + BN1' * ELA * BN1 * detJ * gaussWeight(k) + MAT;
        % Strain energy
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + (a * (trace(CN1) - 3) + b * (trace(GN1) - 3) - d * log(sqrt(I3N1)) + c / 2 * (sqrt(I3N1) - 1)^2) * detJ * gaussWeight(k);
    else
        % stress at gausspoint
        PN1 = FN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FN1) * PN1 * FN1';
        array = postStressComputation(array, N, k, gaussWeight, detJ, detJN1, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = RX;
    kData{1, 1} = KXX;
end
end

function D = secDiffOperator(A)
D = zeros(6, 6);
D(1, 2) = A(3, 3);
D(1, 3) = A(2, 2);
D(1, 5) = -0.5 * (A(2, 3) + A(3, 2));
D(2, 1) = A(3, 3);
D(2, 3) = A(1, 1);
D(2, 6) = -0.5 * (A(1, 3) + A(3, 1));
D(3, 1) = A(2, 2);
D(3, 2) = A(1, 1);
D(3, 4) = -0.5 * (A(1, 2) + A(2, 1));
D(4, 3) = -0.5 * (A(1, 2) + A(2, 1));
D(4, 4) = -0.5 * A(3, 3);
D(4, 5) = 0.25 * (A(1, 3) + A(3, 1));
D(4, 6) = 0.25 * (A(2, 3) + A(3, 2));
D(5, 1) = -0.5 * (A(2, 3) + A(3, 2));
D(5, 4) = 0.25 * (A(1, 3) + A(3, 1));
D(5, 5) = -0.5 * A(1, 1);
D(5, 6) = 0.25 * (A(1, 2) + A(2, 1));
D(6, 2) = -0.5 * (A(1, 3) + A(3, 1));
D(6, 4) = 0.25 * (A(2, 3) + A(3, 2));
D(6, 5) = 0.25 * (A(1, 2) + A(2, 1));
D(6, 6) = -0.5 * A(2, 2);
end