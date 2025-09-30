function [rData, kData, elementEnergy, array] = displacementNeoHookeGENERICMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

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
lambda  = materialObject.lambda;
mu      = materialObject.mu;
bulk    = lambda + 2 / 3 * mu; % bulk modulus
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
%     CxN1inv     = CxN1  \ eye(dimension);
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
    
    if ~computePostData
    %% residual, energy and tangent
        
        % free energy (psi(C, theta) = psi0(theta) + psi1(C) - (theta - thetaR) * psi2(C))
            psi0N        = kappa * (thetaN   - thetaR - thetaN   * log(thetaN   / thetaR));
            psi0N1       = kappa * (thetaN1  - thetaR - thetaN1  * log(thetaN1  / thetaR));
%             psi0N05      = kappa * (thetaN05 - thetaR - thetaN05 * log(thetaN05 / thetaR));
            psiDevN      = mu / 2 * (trace(CxN  ) - 3 - 2 * log(JN  ) - 2 / 3 * (JN   - 1)^2);
            psiDevN1     = mu / 2 * (trace(CxN1 ) - 3 - 2 * log(JN1 ) - 2 / 3 * (JN1  - 1)^2);
%             psiDevN05    = mu / 2 * (trace(CxN05) - 3 - 2 * log(JN05) - 2 / 3 * (JN05 - 1)^2);
            psiVolN      = bulk / 4 * (log(JN  )^2 + (JN   - 1)^2);
            psiVolN1     = bulk / 4 * (log(JN1 )^2 + (JN1  - 1)^2);
%             psiVolN05    = bulk / 4 * (log(JN05)^2 + (JN05 - 1)^2);
            psi1N        = psiDevN   + psiVolN;
            psi1N1       = psiDevN1  + psiVolN1;
%             psi1N05      = psiDevN05 + psiVolN05;
            DpsiVolDJN   = bulk / 2 * (log(JN  ) / JN   + JN   - 1);
            DpsiVolDJN1  = bulk / 2 * (log(JN1 ) / JN1  + JN1  - 1);
            DpsiVolDJN05 = bulk / 2 * (log(JN05) / JN05 + JN05 - 1);
            psi2N        = 3 * beta * DpsiVolDJN;
            psi2N1       = 3 * beta * DpsiVolDJN1;
%             psi2N05      = 3 * beta * DpsiVolDJN05;
            
        % first derivatives
            Dpsi0DthetaN        = -kappa * log(thetaN   / thetaR);
            Dpsi0DthetaN1       = -kappa * log(thetaN1  / thetaR);
%             Dpsi0DthetaN05      = -kappa * log(thetaN05 / thetaR);
%             DDpsi0DthetathetaN   = -kappa / thetaN;
%             DDpsi0DthetathetaN1  = -kappa / thetaN1;
%             DDpsi0DthetathetaN05 = -kappa / thetaN05;
%             DDpsiVolDJJN        = 1 / 2 * bulk * 1 / (JN^2  ) * (1 - log(JN  ) + JN^2  );
%             DDpsiVolDJJN1       = 1 / 2 * bulk * 1 / (JN1^2 ) * (1 - log(JN1 ) + JN1^2 );
            DDpsiVolDJJN05      = 1 / 2 * bulk * 1 / (JN05^2) * (1 - log(JN05) + JN05^2);
%             DpsiDevDCxN         = mu / 2 * (eye(dimension) - CxNinv   * (1 + 2 / 3 * (JN^2   - JN  )));
%             DpsiDevDCxN1        = mu / 2 * (eye(dimension) - CxN1inv  * (1 + 2 / 3 * (JN1^2  - JN1 )));
            DpsiDevDCxN05       = mu / 2 * (eye(dimension) - CxN05inv * (1 + 2 / 3 * (JN05^2 - JN05)));
%             Dpsi2DJN            = 3 * beta * DDpsiVolDJJN;
%             Dpsi2DJN1           = 3 * beta * DDpsiVolDJJN1;
            Dpsi2DJN05          = 3 * beta * DDpsiVolDJJN05;
            
        % internal energy and entropy
            internalEnergyN   = psi0N   - thetaN   * Dpsi0DthetaN   + psi1N   + thetaR * psi2N;
            internalEnergyN1  = psi0N1  - thetaN1  * Dpsi0DthetaN1  + psi1N1  + thetaR * psi2N1;
%             internalEnergyN05 = psi0N05 - thetaN05 * Dpsi0DthetaN05 + psi1N05 + thetaR * psi2N05;
%             entropyN          = -Dpsi0DthetaN   + psi2N;
%             entropyN1         = -Dpsi0DthetaN1  + psi2N1;
%             entropyN05        = -Dpsi0DthetaN05 + psi2N05;
            
%             DentropyDthetaN         = kappa / thetaN;
%             DentropyDthetaN1        = kappa / thetaN1;
            DentropyDthetaN05       = kappa / thetaN05;
%             DentropyDFxN            = 3 * beta * DDpsiVolDJJN   * JN   * FxNinv.';
%             DentropyDFxN1           = 3 * beta * DDpsiVolDJJN1  * JN1  * FxN1inv.';
            DentropyDFxN05          = 3 * beta * DDpsiVolDJJN05 * JN05 * FxN05inv.';
%             DentropyDCxN            = Dpsi2DJN   * 1 / 2 * JN   * CxNinv;
%             DentropyDCxN1           = Dpsi2DJN1  * 1 / 2 * JN1  * CxN1inv;
            DentropyDCxN05          = Dpsi2DJN05 * 1 / 2 * JN05 * CxN05inv;
%             DinternalEnergyDCxN     = DpsiDevDCxN   + (DpsiVolDJN   + thetaR * Dpsi2DJN  ) * 1 / 2 * JN   * CxNinv;
%             DinternalEnergyDCxN1    = DpsiDevDCxN1  + (DpsiVolDJN1  + thetaR * Dpsi2DJN1 ) * 1 / 2 * JN1  * CxN1inv;
            DinternalEnergyDCxN05   = DpsiDevDCxN05 + (DpsiVolDJN05 + thetaR * Dpsi2DJN05) * 1 / 2 * JN05 * CxN05inv;
            
            DentropyDFxN05V = matrixToVoigtNonSymmetric(DentropyDFxN05);
            
        % 2nd Piola Kirchhoff stress tensor        
        SN05    = 2 * (DinternalEnergyDCxN05 - thetaN05 * DentropyDCxN05);
        SN05V   = matrixToVoigt(SN05, 'stress');
        
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
        epsilon = B * uN05(:);
        sigma_V = C * epsilon;
        sigma = voigtToMatrix(sigma_V, 'stress');
        stressTensor.Cauchy = sigma;
        array = postStressComputation(array, N, k, gaussWeight, detJStruct, stressTensor, setupObject, dimension);
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

function [f] = auxilaryFunctionDW(J, lambda)
f = lambda / 2 * (log(J) + J^2 - J);
end
