function [rData, kData, elementEnergy, array] = displacementMooneyRivlinGENERICMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
% 
% Description:
%   material model: Mooney Rivlin in a, b, c1, c2, d1 and d2 plus thermal evolution
%   equation obtained from GENERIC (without L2-smoothing)
%   description in S, C and theta
%   Integrator: Midpoint Rule without L2-smoothing
%
% CREATOR(S)
% 15.01.2025 Moritz Hille, Tim Pleßke

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
dN = shapeFunctionObject.dN_xi_k_I;

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
qN1         = dofs.edN1;
edR         = qR(edof, 1:dimension)';
edN         = qN(edof, 1:dimension)';
edN1        = qN1;
edN05       = 1 / 2 * (edN1 + edN);
vN05        = (edN1 - edN) / DT;
tempN1      = dofs.thetaN1;
tempN       = qN(edof,dimension+1).';
tempN05     = 1 / 2 * (tempN + tempN1);

%% material data

a       = materialObject.a;
b       = materialObject.b;
c1      = materialObject.c1;
c2      = materialObject.c2;
d1      = materialObject.d1;
d2      = materialObject.d2;

kappa   = materialObject.kappa;
beta    = materialObject.beta;
thetaR  = materialObject.thetaR;
k0      = materialObject.k0;


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

%% Run through all Gauss points
for k = 1:numberOfGausspoints

    % Jacobian matrix and determinant
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR, edN, edN1, dN, k, setupObject);

    % deformation gradient
    gradvN05    = vN05 * dN_X_I.';
    gradvN05V   = matrixToVoigtNonSymmetric(gradvN05);

    FxN1        = edN1  * dN_X_I.';
    FxN         = edN   * dN_X_I.';
    FxN05       = edN05 * dN_X_I.';
    FxN05inv    = FxN05 \ eye(dimension);
    CxN1        = FxN1.'  * FxN1;
    CxN         = FxN.'   * FxN;
    CxN05       = FxN05.' * FxN05;
    CxN1inv     = CxN1 \ eye(dimension);
    CxNinv      = CxN \ eye(dimension);
    CxN05inv    = CxN05 \ eye(dimension);


    % B-matrix (midpoint)
    BN05 = BMatrix(dN_X_I, FxN05);

    % auxiliary variables
    JN05 = sqrt(det(CxN05));
    JN1  = sqrt(det(CxN1 ));
    JN   = sqrt(det(CxN  ));

    % temperature
    thetaN1         = N(k, :) * tempN1.';
    thetaN          = N(k, :) * tempN.';
    thetaN05        = 1 / 2 * (thetaN + thetaN1);
    gradthetaN05    = tempN05 * dN_X_I.';

        %% residual, energy and tangent

        % free energy (psi(C, G, J, theta) = psi0(theta) + psi1(C, G, J) + (theta - thetaR) * psi2(C))
        psi0N        = kappa * (thetaN   - thetaR - thetaN   * log(thetaN   / thetaR));
        psi0N1       = kappa * (thetaN1  - thetaR - thetaN1  * log(thetaN1  / thetaR));
        psi0N05      = kappa * (thetaN05 - thetaR - thetaN05 * log(thetaN05 / thetaR));

        GxN          = JN^2*CxNinv;
        GxN1         = JN1^2*CxN1inv;
        GxN05        = JN05^2*CxN05inv;
        
        psi1aN       = a * (trace(CxN) - 3) + b * (trace(GxN) - 3);
        psi1aN1      = a * (trace(CxN1) - 3) + b * (trace(GxN1) - 3);
        psi1aN05     = a * (trace(CxN05) - 3) + b * (trace(GxN05) - 3);

        psi1bN       = c1/2 * (JN - 1)^2 - d1 * log(JN);
        psi1bN1      = c1/2 * (JN1 - 1)^2 - d1 * log(JN1);
        psi1bN05     = c1/2 * (JN05 - 1)^2 - d1 * log(JN05);

        psi1N        = psi1aN   + psi1bN;
        psi1N1       = psi1aN1  + psi1bN1;
        psi1N05      = psi1aN05 + psi1bN05;
        
        psi2N        = 3 * beta * (c2 * (JN-1) - d2*JN^(-1));
        psi2N1       = 3 * beta * (c2 * (JN1-1) - d2*JN1^(-1));
        psi2N05      = 3 * beta * (c2 * (JN05-1) - d2*JN05^(-1));

        % first derivatives
        Dpsi0DthetaN        = -kappa * log(thetaN   / thetaR);
        Dpsi0DthetaN1       = -kappa * log(thetaN1  / thetaR);
        Dpsi0DthetaN05      = -kappa * log(thetaN05 / thetaR);

        % DDpsi0DthetathetaN   = -kappa / thetaN;
        % DDpsi0DthetathetaN1  = -kappa / thetaN1;
        % DDpsi0DthetathetaN05 = -kappa / thetaN05;

        Dpsi2DJN            = 3 * beta * (c2 + d2*JN^(-2));
        Dpsi2DJN1           = 3 * beta * (c2 + d2*JN1^(-2));
        Dpsi2DJN05          = 3 * beta * (c2 + d2*JN05^(-2));

        Dpsi1aDCxN = (a + b*trace(CxN))* eye(dimension) - b * CxN;
        Dpsi1aDCxN1 = (a + b*trace(CxN1))* eye(dimension) - b * CxN1;
        Dpsi1aDCxN05 = (a + b*trace(CxN05))* eye(dimension) - b * CxN05;

        Dpsi1bDCxN  = (c1*(JN -1) - d1*JN^(-1))*1/2*JN*CxNinv;
        Dpsi1bDCxN1  = (c1*(JN1 -1) - d1/JN1^(-1))*1/2*JN1*CxN1inv;
        Dpsi1bDCxN05  = (c1*(JN05 -1) - d1*JN05^(-1))*1/2*JN05*CxN05inv;

        Dpsi1DCxN    = Dpsi1aDCxN + Dpsi1bDCxN ;
        Dpsi1DCxN1    = Dpsi1aDCxN1 + Dpsi1bDCxN1 ;
        Dpsi1DCxN05    = Dpsi1aDCxN05 + Dpsi1bDCxN05 ;
        




        % internal energy and entropy
        internalEnergyN         = psi0N   - thetaN   * Dpsi0DthetaN   + psi1N   - thetaR * psi2N;
        internalEnergyN1        = psi0N1  - thetaN1  * Dpsi0DthetaN1  + psi1N1  - thetaR * psi2N1;
        % internalEnergyN05       = psi0N05 - thetaN05 * Dpsi0DthetaN05 + psi1N05 - thetaR * psi2N05;

        % entropyN                = -Dpsi0DthetaN   - psi2N;
        % entropyN1               = -Dpsi0DthetaN1  - psi2N1;
        % entropyN05              = -Dpsi0DthetaN05 - psi2N05;

        % DentropyDthetaN         = kappa / thetaN;
        % DentropyDthetaN1        = kappa / thetaN1;
        DentropyDthetaN05       = kappa / thetaN05;

        % DentropyDFxN            = -Dpsi2DJN   * JN   * FxNinv.';
        % DentropyDFxN1           = -Dpsi2DJN1  * JN1  * FxN1inv.';
        DentropyDFxN05          = -Dpsi2DJN05 * JN05 * FxN05inv.';

%         DentropyDCxN            = -Dpsi2DJN   * 1 / 2 * JN   * CxNinv;
        DentropyDCxN1           = -Dpsi2DJN1  * 1 / 2 * JN1  * CxN1inv;
        DentropyDCxN05          = -Dpsi2DJN05 * 1 / 2 * JN05 * CxN05inv;

        % DinternalEnergyDCxN     = Dpsi1DCxN   - thetaR * Dpsi2DJN  * 1 / 2 * JN   * CxNinv;
        DinternalEnergyDCxN1    = Dpsi1DCxN1  - thetaR * Dpsi2DJN1  * 1 / 2 * JN1  * CxN1inv;
        DinternalEnergyDCxN05   = Dpsi1DCxN05 - thetaR * Dpsi2DJN05 * 1 / 2 * JN05 * CxN05inv;

        DentropyDFxN05V = matrixToVoigtNonSymmetric(DentropyDFxN05);

        % 2nd Piola Kirchhoff stress tensor
        SN05   = 2 * (DinternalEnergyDCxN05 - thetaN05 * DentropyDCxN05);
        SN1    = 2 * (DinternalEnergyDCxN1  - thetaN1  * DentropyDCxN1 );        
        
        SN05V   = matrixToVoigt(SN05, 'stress');
        
    if ~computePostData

        % heat flux vector
        QN05 = - k0 * JN05 * CxN05inv * gradthetaN05.';

        % residual
        RX = RX + BN05.' * SN05V * detJ * gaussWeight(k);
        RT = RT + (N(k,:).' * 1 / DT * (thetaN1 - thetaN) +  N(k,:).' * 1 / DentropyDthetaN05 * (DentropyDFxN05V.' * gradvN05V) - 1 / kappa * dN_X_I.' * QN05) * detJ * gaussWeight(k);

        % tangent parts

        % energy and energy difference
        elementEnergy.internalEnergy = elementEnergy.internalEnergy + internalEnergyN1 * detJ * gaussWeight(k);
        elementEnergy.deltaU = elementEnergy.deltaU + (internalEnergyN1 - internalEnergyN) * detJ * gaussWeight(k);

        % tangent



    else
        % stress at gausspoint
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N, k, gaussWeight, detJ, stressTensor, setupObject, dimension);

    end
end
% if setupObject.newton.step(setupObject.timeStep) == 2
%    keyboard;
% end
if ~computePostData
    rData{1} = RX;
    rData{2} = RT;
    kData{1, 1} = KXX;
    kData{1, 2} = KXT;
    kData{2, 1} = KTX;
    kData{2, 2} = KTT;
end
end

function Av = matrixToVoigtNonSymmetric(A)
%% Converts nonsymmetric matrix into voigt notation
Av = [A(1,1); A(2,2); A(3,3); A(1,2); A(2,3); A(1,3); A(2,1); A(3,2); A(3,1)];
end




