function [rData, kData, elementEnergy, array] = displacementNeoHookeGENERICEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
%
% Description:
%   material model: Neo Hooke in lambda and mu plus thermal evolution equation obtained from GENERIC
%   description in S, C and theta
%   Integrator: Midpoint Rule
%
% 16.05.2024 Moritz Hille

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;

% element degree of freedom tables and more
edof = meshObject.edof(e, :);
dimension = obj.dimension;
DT = setupObject.timeStepSize;

% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

% Voigt maps
% map.Voigt = [1, 5, 9, 4, 8, 7]';
% map.VoigtInv = [1, 4, 6; 4, 2, 5; 6, 5, 3];
% map.VoigtFull = [1, 5, 9, 4, 8, 7, 2, 6, 3]';
% map.VoigtFullInv = [1, 4, 6; 7, 2, 5; 9, 8, 3];
% map.Isym = [eye(3), zeros(3); zeros(3), 2 * eye(3)];
% map.Isyminv = [eye(3), zeros(3); zeros(3), 1 / 2 * eye(3)];
% mapVoigtObject = obj.mapVoigtObject;

%% extract dofs for element
% nodal dofs
qR          = obj.qR;
qN          = obj.qN;
qN1         = obj.qN1;
edR         = qR(edof, 1:dimension)';
edN         = qN(edof, 1:dimension)';
edN1        = qN1(edof, 1:dimension)';
velocityN   = obj.vN(edof, 1:dimension)';
tempN1      = qN1(edof,dimension+1).';
tempN       = qN(edof,dimension+1).';

%% material data
rho     = materialObject.rho;
lambda  = materialObject.lambda;
mu      = materialObject.mu;
kappa   = materialObject.kappa;
beta    = materialObject.beta;
thetaR  = materialObject.thetaR;
k0      = materialObject.k0;

a      = materialObject.a;
b      = materialObject.b;
c      = materialObject.c;
d      = materialObject.d;


%% initialize energy, residual and tangent (global)

% initialize internal energy and difference of internal energy
elementEnergy.internalEnergy = 0;
elementEnergy.deltaU = 0;
elementEnergy.helmholtzEnergy = 0;

% Initialize sub-residual vector
RX = rData{1};
RT = rData{2};

% Initialize sub-tangent matrices (X = displacement part, T = temperature part)
KXX = kData{1, 1};
KXT = kData{1, 2};
KTX = kData{2, 1};
KTT = kData{2, 2};
%     dimR = numel(rData{1});

% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

%% Run through all Gauss points
for k = 1:numberOfGausspoints
    % Jacobian matrix and determinant
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    FxN        = edN * dN_X_I';
    FxN1        = edN1 * dN_X_I';
    % FxN1inv    = FxN1 \ eye(dimension);
    CxN1        = FxN1.' * FxN1;
    CxN1inv    = CxN1 \ eye(dimension);

    % B-matrix (midpoint)
    BN1    = BMatrix(dN_X_I, FxN1);

    % auxiliary variables
    JN1     = det(FxN1);

    if ~computePostData


        %% residual, energy and tangent
        if 0
            % Psi = mu/2 * (trace(C) - 3 - 2*log(J) - 2/3 *(J - 1)^2) + 1/4 * (lambda * 2/3 * mu) * (log(J)^2 + (J - 1)^2);
            dPsi_J = mu/2*(-2/JN1-4/3*(JN1-1)) + 1/4*(lambda+2/3*mu)*(2*log(JN1)*1/JN1+2*(JN1-1));
            dJ_C = 1/2*JN1*CxN1inv;
            dPsi_C = mu/2*eye(3);
        else
            % Psi = a*(trace(CxN1)-3) + b*(trace(GxN1)-3) - d*log(JN1) + c/2*(JN1-1)^2;
            dPsi_J = d/JN1 + c*(JN1-1);
            dJ_C = 1/2*JN1*CxN1inv;
            dPsi_C = a*eye(3);
        end
        SN1 = 2*(dPsi_C + dPsi_J*dJ_C);
        SN1V   = matrixToVoigt(SN1, 'stress');
        RX = RX + BN1.' * SN1V * detJ * gaussWeight(k);

    else
        % stress at gausspoint
        epsilon = B * uN05(:);
        sigma_V = C * epsilon;
        sigma = voigtToMatrix(sigma_V, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array, N, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    rData{1} = RX;
    rData{2} = RT;
    kData{1, 1} = KXX;
    kData{1, 2} = KXT;
    kData{2, 1} = KTX;
    kData{2, 2} = KTT;
end
end
