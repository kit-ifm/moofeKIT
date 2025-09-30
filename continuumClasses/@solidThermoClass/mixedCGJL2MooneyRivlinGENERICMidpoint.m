function [rData, kData, elementEnergy, array] = mixedCGJL2MooneyRivlinGENERICMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
%
% Description:
%   material model: Mooney Rivlin in C, G and J plus thermal evolution equation obtained from GENERIC
%   description in S, and mixed variables C, G and J aswell as theta with an additional mixed projection method for partial
%   derivatives of internal energy and entropy with respect to theta
%   Integrator: Midpoint Rule with L2-smoothing
%
% 25.11.2024 Moritz Hille

%% setup
    % load objects
    shapeFunctionObject = obj.shapeFunctionObject;
    materialObject = obj.materialObject;
    mixedFEObject = obj.mixedFEObject;
    meshObject = obj.meshObject;
    
    % element degree of freedom tables and more
    edof = meshObject.edof(e, :);
    dimension = obj.dimension;
    DT = setupObject.timeStepSize;

    % gauss integration and shape functions
    gaussWeight = shapeFunctionObject.gaussWeight;
    numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
    N  = shapeFunctionObject.N_k_I;
    dN = shapeFunctionObject.dN_xi_k_I;
    M  = mixedFEObject.shapeFunctionObject.N_k_I;
    dM = mixedFEObject.shapeFunctionObject.dN_xi_k_I;
    numberOfInternalNodes = size(M, 2);

%% extract dofs for element
    % nodal dofs
    qR          = obj.qR;
    qN          = obj.qN;
    qN1         = dofs.edN1;
    edR         = qR(edof, 1:dimension)';
    edN         = qN(edof, 1:dimension)';
    edN1        = qN1;
    edN05       = 1 / 2 * (edN1 + edN);
%     vN05        = (edN1 - edN) / DT;
    tempN1      = dofs.thetaN1;
    tempN       = qN(edof,dimension+1).';
    tempN05     = 1 / 2 * (tempN + tempN1);
    edAlphaN    = obj.mixedFEObject.qN(e, :);
    edAlphaN1   = dofs.edAlphaN1;
    edAlphaN05  = 1 / 2 * (edAlphaN + edAlphaN1);
    
    % mixed dofs
%     mixedEN     = edAlphaN(1:numberOfInternalNodes);
%     mixedEtaN   = edAlphaN(numberOfInternalNodes+1:2*numberOfInternalNodes);
    mixedCNV    = edAlphaN(2*numberOfInternalNodes +1:8*numberOfInternalNodes);
    mixedGNV    = edAlphaN(8*numberOfInternalNodes +1:14*numberOfInternalNodes);
    mixedJN     = edAlphaN(14*numberOfInternalNodes +1:15*numberOfInternalNodes);

%     mixedEN1    = edAlphaN1(1:numberOfInternalNodes);
%     mixedEtaN1  = edAlphaN1(numberOfInternalNodes+1:2*numberOfInternalNodes);
    mixedCN1V   = edAlphaN1(2*numberOfInternalNodes +1:8*numberOfInternalNodes);
    mixedGN1V   = edAlphaN1(8*numberOfInternalNodes +1:14*numberOfInternalNodes);
    mixedJN1    = edAlphaN1(14*numberOfInternalNodes +1:15*numberOfInternalNodes);

    mixedEN05   = edAlphaN05(1:numberOfInternalNodes);
    mixedEtaN05 = edAlphaN05(numberOfInternalNodes+1:2*numberOfInternalNodes);
    mixedCN05V  = edAlphaN05(2*numberOfInternalNodes +1:8*numberOfInternalNodes);
    mixedGN05V  = edAlphaN05(8*numberOfInternalNodes +1:14*numberOfInternalNodes);
    mixedJN05   = edAlphaN05(14*numberOfInternalNodes +1:15*numberOfInternalNodes);
    

%% material data
matA  = materialObject.a;
matB  = materialObject.b;
matC1 = materialObject.c1;
matD1 = materialObject.d1;
matC2 = materialObject.c2;
matD2 = materialObject.d2;

% lambda  = materialObject.lambda;
% mu      = materialObject.mu;
% bulk    = lambda + 2 / 3 * mu; % bulk modulus
kappa   = materialObject.kappa;
beta    = materialObject.beta;
thetaR  = materialObject.thetaR;
k0      = materialObject.k0;


%% initialize energy, residual and tangent (global)
    
    % initialize internal energy and difference of internal energy
    elementEnergy.internalEnergy = 0;
    elementEnergy.deltaU = 0;
    elementEnergy.helmholtzEnergy = 0;
    elementEnergy.entropy = 0;
    elementEnergy.deltaS = 0;
    
    % Initialize sub-residual vector
    RX   = rData{1};
    RT   = rData{2};
    RE   = rData{3};
    REta = rData{4};
    RC   = rData{5};
    RG   = rData{6};
    RJ   = rData{7};
    
    % Initialize sub-tangent matrices (X = displacement part, T = temperature part)
    KXX   = kData{1, 1}; KXT   = kData{1, 2}; KXE   = kData{1, 3}; KXEta   = kData{1, 4}; KXC     = kData{1, 5}; KXG     = kData{1, 6}; KXJ     = kData{1, 7};
    KTX   = kData{2, 1}; KTT   = kData{2, 2}; KTE   = kData{2, 3}; KTEta   = kData{2, 4}; KTC     = kData{2, 5}; KTG     = kData{2, 6}; KTJ     = kData{2, 7};
    KEX   = kData{3, 1}; KET   = kData{3, 2}; KEE   = kData{3, 3}; KEEta   = kData{3, 4}; KEC     = kData{3, 5}; KEG     = kData{3, 6}; KEJ     = kData{3, 7};
    KEtaX = kData{4, 1}; KEtaT = kData{4, 2}; KEtaE = kData{4, 3}; KEtaEta = kData{4, 4}; KEtaC   = kData{4, 5}; KEtaG   = kData{4, 6}; KEtaJ   = kData{4, 7};
    KCX   = kData{5, 1}; KCT   = kData{5, 2}; KCE   = kData{5, 3}; KCEta   = kData{5, 4}; KCC     = kData{5, 5}; KCG     = kData{5, 6}; KCJ     = kData{5, 7};
    KGX   = kData{6, 1}; KGT   = kData{6, 2}; KGE   = kData{6, 3}; KGEta   = kData{6, 4}; KGC     = kData{6, 5}; KGG     = kData{6, 6}; KGJ     = kData{6, 7};
    KJX   = kData{7, 1}; KJT   = kData{7, 2}; KJE   = kData{7, 3}; KJEta   = kData{7, 4}; KJC     = kData{7, 5}; KJG     = kData{7, 6}; KJJ     = kData{7, 7};
    
%% Run through all Gauss points
for k = 1:numberOfGausspoints
    
    % Jacobian matrix and determinant
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR, edN, edN1, dN, k, setupObject);
    [~, ~, dM_X_I, ~, ~, ~] = computeAllJacobian(edR, edN, edN1, dM, k, setupObject);

    % deformation gradient
    FxN1            = edN1  * dN_X_I.';
    FxN             = edN   * dN_X_I.';
    FxN05           = edN05 * dN_X_I.';
%     FxNV            = matrixToVoigtUnsymmetric(FxN);
%     FxN1V           = matrixToVoigtUnsymmetric(FxN1);
%     FxN05V          = matrixToVoigtUnsymmetric(FxN05);
%     FxN05inv        = FxN05 \ eye(dimension);
%     FxN05invV       = matrixToVoigtUnsymmetric(FxN05inv);
%     FxN05invtranspV = matrixToVoigtUnsymmetric(FxN05inv.');
 
%     CxN1            = FxN1.'  * FxN1;
%     CxN             = FxN.'   * FxN;
%     CxN05           = FxN05.' * FxN05;
%     CxN05inv        = CxN05 \ eye(dimension);
%     CxN05invV       = matrixToVoigt(CxN05inv, 'stress');
        
    % B-matrix (midpoint)
    BN05 = BMatrix(dN_X_I, FxN05);
    BN1  = BMatrix(dN_X_I, FxN1);
    
    % temperature
    thetaN1         = N(k, :) * tempN1.';
    thetaN          = N(k, :) * tempN.';
    thetaN05        = 1 / 2 * (thetaN + thetaN1);
    gradthetaN05    = tempN05 * dN_X_I.';
    
    % mixed variables
    EN05     = M(k,:) * mixedEN05.';
    EtaN05   = M(k,:) * mixedEtaN05.';
    
    CNV   = reshape(mixedCNV, 6, []) * M(k, :)';
    CN    = voigtToMatrix(CNV, 'strain');
    GNV   = reshape(mixedGNV, 6, []) * M(k, :)';
    GN    = voigtToMatrix(GNV, 'strain');
    JN    = M(k, :) * mixedJN.';
    CN1V  = reshape(mixedCN1V, 6, []) * M(k, :)';
    CN1   = voigtToMatrix(CN1V, 'strain');
    GN1V  = reshape(mixedGN1V, 6, []) * M(k, :)';
    GN1   = voigtToMatrix(GN1V, 'strain');
    JN1   = M(k, :) * mixedJN1.';
    CN05V = reshape(mixedCN05V, 6, []) * M(k, :)';
    CN05  = voigtToMatrix(CN05V, 'strain');
    GN05V = reshape(mixedGN05V, 6, []) * M(k, :)';
    GN05  = voigtToMatrix(GN05V, 'strain');
    JN05  = M(k, :) * mixedJN05.';
    
    I = eye(dimension);
    CN05inv  = CN05\I;
    CN05invV = matrixToVoigt(CN05inv, 'stress');
    
    gradEN05 = mixedEN05 * dM_X_I.';
    
    %% residual, energy and tangent
        
        % free energy (W(C,G,J,theta) = Psi0(theta) + Psi1a(C,G) + Psi1b(J) + (theta - thetaR) * Psi2(J))
            psi0N    = kappa * (thetaN   - thetaR - thetaN   * log(thetaN   / thetaR));
            psi0N1   = kappa * (thetaN1  - thetaR - thetaN1  * log(thetaN1  / thetaR));
%             psi0N05  = kappa * (thetaN05 - thetaR - thetaN05 * log(thetaN05 / thetaR));
            psi1aN   = matA * (trace(CN  ) - 3) + matB * (trace(GN  ) - 3);
            psi1aN1  = matA * (trace(CN1 ) - 3) + matB * (trace(GN1 ) - 3);
%             psi1aN05 = matA * (trace(CN05) - 3) + matB * (trace(GN05) - 3);
            psi1bN   = matC1 / 2 * (JN   - 1)^2 - matD1 * log(JN  );
            psi1bN1  = matC1 / 2 * (JN1  - 1)^2 - matD1 * log(JN1 );
%             psi1bN05 = matC1 / 2 * (JN05 - 1)^2 - matD1 * log(JN05);
            psi1N    = psi1aN   + psi1bN;
            psi1N1   = psi1aN1  + psi1bN1;
%             psi1N05  = psi1aN05 + psi1bN05;
            psi2N    = 3 * beta * (matC2 * (JN   - 1) - matD2 * 1 / JN  );
            psi2N1   = 3 * beta * (matC2 * (JN1  - 1) - matD2 * 1 / JN1 );
%             psi2N05  = 3 * beta * (matC2 * (JN05 - 1) - matD2 * 1 / JN05);
            
        % derivatives 
            Dpsi0DthetaN          = -kappa * log(thetaN   / thetaR);
            Dpsi0DthetaN1         = -kappa * log(thetaN1  / thetaR);
%             Dpsi0DthetaN05        = -kappa * log(thetaN05 / thetaR);
%             Dpsi1aDCN             = matA * I;
            Dpsi1aDCN1            = matA * I;
            Dpsi1aDCN05           = matA * I;
%             Dpsi1aDGN             = matB * I;
            Dpsi1aDGN1            = matB * I;
            Dpsi1aDGN05           = matB * I;
%             Dpsi1bDJN             = matC1 * (JN   - 1) - matD1 * 1 / JN;
            Dpsi1bDJN1            = matC1 * (JN1  - 1) - matD1 * 1 / JN1;
            Dpsi1bDJN05           = matC1 * (JN05 - 1) - matD1 * 1 / JN05;
%             Dpsi2DJN              = 3 * beta * (matC2 + matD2 * 1 / (JN  ^2));
            Dpsi2DJN1             = 3 * beta * (matC2 + matD2 * 1 / (JN1 ^2));
            Dpsi2DJN05            = 3 * beta * (matC2 + matD2 * 1 / (JN05^2));
%             DDpsi0DthetaDthetaN   = - kappa / thetaN; 
%             DDpsi0DthetaDthetaN1  = - kappa / thetaN1;
            DDpsi0DthetaDthetaN05 = - kappa / thetaN05;
%             DDpsi1bDJDJN          = matC1 + matD2 * 1 / (JN  ^2);
%             DDpsi1bDJDJN1         = matC1 + matD2 * 1 / (JN1 ^2);
            DDpsi1bDJDJN05        = matC1 + matD1 * 1 / (JN05^2);
%             DDpsi2DJDJN           = - 6 * beta * matD2 * 1 / (JN  ^3);
%             DDpsi2DJDJN1          = - 6 * beta * matD2 * 1 / (JN1 ^3);
            DDpsi2DJDJN05         = - 6 * beta * matD2 * 1 / (JN05^3);
          
        % internal energy and entropy
            internalEnergyN   = psi0N   - thetaN   * Dpsi0DthetaN   + psi1N   - thetaR * psi2N;
            internalEnergyN1  = psi0N1  - thetaN1  * Dpsi0DthetaN1  + psi1N1  - thetaR * psi2N1;
%             internalEnergyN05 = psi0N05 - thetaN05 * Dpsi0DthetaN05 + psi1N05 - thetaR * psi2N05;

            DinternalEnergyDthetaN05    = - thetaN05 * DDpsi0DthetaDthetaN05;
            DinternalEnergyDCN1         = Dpsi1aDCN1;
            DinternalEnergyDCN05        = Dpsi1aDCN05;
            DinternalEnergyDGN1         = Dpsi1aDGN1;
            DinternalEnergyDGN05        = Dpsi1aDGN05;
            DinternalEnergyDJN1         = Dpsi1bDJN1  - thetaR * Dpsi2DJN1;
            DinternalEnergyDJN05        = Dpsi1bDJN05 - thetaR * Dpsi2DJN05;
            
            DDinternalEnergyDJDJN05 = DDpsi1bDJDJN05 - thetaR * DDpsi2DJDJN05; 
            
            entropyN            = -Dpsi0DthetaN - psi2N;
            entropyN1           = -Dpsi0DthetaN1 - psi2N1;
%             entropyN05          = - DWDthetaN05;

%             DentropyDthetaN     = -DDpsi0DthetaDthetaN;
%             DentropyDthetaN1    = -DDpsi0DthetaDthetaN1;
            DentropyDthetaN05   = -DDpsi0DthetaDthetaN05;
%             DentropyDJN         = - Dpsi2DJN;
            DentropyDJN1        = - Dpsi2DJN1;
            DentropyDJN05       = - Dpsi2DJN05;

            DDentropyDJDJN05    = -DDpsi2DJDJN05;

             

            
%             DentropyDthetaN05       = kappa / thetaN05;
%             DentropyDFxN05          = 3 * beta * DDpsiVolDJJN05 * JN05 * FxN05inv.';
%             DentropyDCxN05          = Dpsi2DJN05 * 1 / 2 * JN05 * CxN05inv;
%             DinternalEnergyDCxN05   = DpsiDevDCxN05 + (DpsiVolDJN05 + thetaR * Dpsi2DJN05) * 1 / 2 * JN05 * CxN05inv;
%             DentropyDFxN05V         = matrixToVoigtUnsymmetric(DentropyDFxN05);
            
%             DDinternalEnergyDCxCxN05V   = DDpsiDevDCxCxN05V + (DDpsiVolDJJN05 + thetaR * DDpsi2DJJN05) * 1 / 4 * JN05^2 * (CxN05invV * CxN05invV.') + 1 / 2 * (DpsiVolDJN05 + thetaR * Dpsi2DJN05) * JN05 * (1 / 2 * (CxN05invV * CxN05invV.') - sym_voigt(CxN05inv, CxN05inv, 'stress', 'stress'));
            
        % 2nd Piola Kirchhoff stress tensor        
            SN1         = 2 * DinternalEnergyDCN1  + 2 * wedge(DinternalEnergyDGN1, CN1 ) + 1 / JN1  * GN1  * DinternalEnergyDJN1  - thetaN1  * 1 / JN1  * GN1 *  DentropyDJN1;
            SN05        = 2 * DinternalEnergyDCN05 + 2 * wedge(DinternalEnergyDGN05,CN05) + 1 / JN05 * GN05 * DinternalEnergyDJN05 - thetaN05 * 1 / JN05 * GN05 * DentropyDJN05;
            SN05V       = matrixToVoigt(SN05, 'stress');
%             DSDCxN05V    = 2 * (DDinternalEnergyDCxCxN05V - thetaN05 * DDentropyDCxCxN05V);
%             DSDthetaN05  = -2 * DentropyDCxN05;
%             DSDthetaN05V = matrixToVoigt(DSDthetaN05, 'stress');
        
    if ~computePostData
        
        % heat flux vector
            QN05 = - k0 * JN05 * CN05inv * gradthetaN05.';
        
        % residual
            DCx  = 1 / DT * (FxN1.'*FxN1 - FxN.'*FxN);
                    
            RX   = RX   + BN05.' * SN05V * detJ * gaussWeight(k);
            RT   = RT   + (N(k,:).' * 1 / DT * (thetaN1 - thetaN) +  N(k,:).' * (1 / EtaN05) * (DentropyDJN05 * (1 / 2) * (1 / JN05) * GN05V.' * matrixToVoigt(DCx, 'stress')) ...
                            - (1 / EN05 * dN_X_I.' - N(k,:).' / EN05^2 * gradEN05) * QN05) * detJ * gaussWeight(k);
            RE   = RE   + M(k,:).' * (DinternalEnergyDthetaN05 - EN05) * detJ * gaussWeight(k);
            REta = REta + M(k,:).' * (DentropyDthetaN05 - EtaN05) * detJ * gaussWeight(k);
            RC   = RC   + kron(M(k,:).', (1 / DT * (CN1V - CNV) - matrixToVoigt(DCx, 'strain'))) * detJ * gaussWeight(k);
            RG   = RG   + kron(M(k,:).', (1 / DT * (GN1V - GNV) - matrixToVoigt(wedge(CN05,DCx),'strain'))) * detJ * gaussWeight(k);
            RJ   = RJ   + M(k,:).' * (1 / DT * (JN1 - JN) - 1 / 2 * 1 / JN05 * GN05V.' * matrixToVoigt(DCx, 'stress')) * detJ * gaussWeight(k);
            
        % tangent
            IV      = [ 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1/2 0 0; 0 0 0 0 1/2 0; 0 0 0 0 0 1/2]; % Changes a Voigt stress vector into a Voigt strain vector
            KXX_GP     = kron(dN_X_I.' * 1 / 2 * SN05 * dN_X_I, eye(dimension));
            KXT_GP     = -BN05.' * 1 / 2 * 1 / JN05 * DentropyDJN05 * matrixToVoigt(GN05,'stress') * N(k,:);
%             KXC_GP     = BN05.' * DiffSec(DinternalEnergyDGN05) * kron(M(k,:), eye(6));
            KXC_GP    = BN05.' * secDiffOperator(DinternalEnergyDGN05) * IV * IV * kron(M(k,:), eye(6)); %????
            KXG_GP     = BN05.' * 1 / 2 * 1 / JN05 * (DinternalEnergyDJN05 - thetaN05 * DentropyDJN05)* IV * kron(M(k,:), eye(6));
            KXJ_GP     = BN05.' * 1 / 2 * (-1 / (JN05^2) * (DinternalEnergyDJN05 - thetaN05 * DentropyDJN05) + 1 / JN05 * (DDinternalEnergyDJDJN05 - thetaN05 * DDentropyDJDJN05))* matrixToVoigt(GN05,'stress') * M(k,:);
                        
            KTX_GP     = N(k,:).' / EtaN05 * 1 / DT * DentropyDJN05 * 1 / JN05 * matrixToVoigt(GN05,'stress').' * BN1; 
            KTT_GP     = 1 / DT * N(k,:).' * N(k,:) + ((1 / EN05 * dN_X_I.') - (N(k,:).' / (EN05^2) * gradEN05)) * k0 * JN05 * CN05inv * 1 / 2 * dN_X_I;
            KTE_GP     = 1 / 2 * 1 / (EN05^2) * dN_X_I.' * QN05 * M(k,:) - N(k,:).' / (EN05^3) * gradEN05 * QN05 * M(k,:) + 1 / 2 * N(k,:).' / (EN05^2) * QN05.' * dM_X_I;
            KTEta_GP   = -1 / 4 * N(k,:).' * 1 / (EtaN05^2) * DentropyDJN05 * 1 / JN05 * (GN05V.' * matrixToVoigt(DCx,'stress')) * M(k,:); 
            KTC_GP     = -1 / 2 * k0 * JN05 * (1 / EN05 * (BMatrixTemp(dN_X_I,gradthetaN05)).' - N(k,:).' / (EN05^2) * matrixToVoigt(1 / 2 * (gradEN05.' * gradthetaN05 + gradthetaN05.' * gradEN05),'strain').')...
                            * sym_voigt(CN05inv, CN05inv, 'stress', 'stress') * kron(M(k,:), eye(6));
            KTG_GP     = 1 / 4 * N(k,:).' / EtaN05 * DentropyDJN05 * 1 / JN05 * (matrixToVoigt(DCx, 'stress').' * kron(M(k,:), eye(6))); %???
            KTJ_GP     = N(k,:).' / EtaN05 * 1 / 4 * (-1 / (JN05^2) * DentropyDJN05 + 1 / JN05 * DDentropyDJDJN05) * (GN05V.' * matrixToVoigt(DCx, 'stress')) * M(k,:) ...
                            + 1 / 2 * k0 * (1 / EN05 * (BMatrixTemp(dN_X_I,gradthetaN05)).' - N(k,:).' / (EN05^2) * matrixToVoigt(1 / 2 * (gradEN05.' * gradthetaN05 + gradthetaN05.' * gradEN05),'strain').') * CN05invV * M(k,:);
            
            KEE_GP     = -1 / 2 * M(k,:).' * M(k,:);
            
            KEtaT_GP   = -1 / 2 * M(k,:).' * kappa / (thetaN05^2) * N(k,:);
            KEtaEta_GP = -1 / 2 * M(k,:).' * M(k,:);
            
            KCX_GP     = - 2 / DT * kron(M(k,:).', eye(6))* BN1; 
            KCC_GP     = 1 / DT * kron((M(k,:).' * M(k,:)), eye(6));
            
            KGX_GP     = -2 / DT * kron(M(k,:).', eye(6)) * secDiffOperator(CN05) * IV *  BN1; 
            KGC_GP     = -1 / 2 * kron(M(k,:).', eye(6)) * secDiffOperator(DCx) * IV * kron(M(k,:), eye(6));
            KGG_GP     = 1 / DT * kron((M(k,:).' * M(k,:)), eye(6));
            
            KJX_GP     = -M(k,:).' * 1 / DT * 1 / JN05 * (matrixToVoigt(GN05,'stress').' * BN1);
            KJG_GP     = -M(k,:).' * 1 / 4 * 1 / JN05 * (matrixToVoigt(DCx, 'stress').' * kron (M(k,:),eye(6)));
            KJJ_GP     = M(k,:).' * (1 / DT  + 1 / 4 * 1 / (JN05^2) * (GN05V.' * matrixToVoigt(DCx, 'stress'))) * M(k,:); 

            KXX     = KXX     + KXX_GP     * detJ * gaussWeight(k);
            KXT     = KXT     + KXT_GP     * detJ * gaussWeight(k);
            KXC     = KXC     + KXC_GP     * detJ * gaussWeight(k);
            KXG     = KXG     + KXG_GP     * detJ * gaussWeight(k);
            KXJ     = KXJ     + KXJ_GP     * detJ * gaussWeight(k);
            KTX     = KTX     + KTX_GP     * detJ * gaussWeight(k);
            KTT     = KTT     + KTT_GP     * detJ * gaussWeight(k);
            KTE     = KTE     + KTE_GP     * detJ * gaussWeight(k);
            KTEta   = KTEta   + KTEta_GP   * detJ * gaussWeight(k);
            KTC     = KTC     + KTC_GP     * detJ * gaussWeight(k);
            KTG     = KTG     + KTG_GP     * detJ * gaussWeight(k);
            KTJ     = KTJ     + KTJ_GP     * detJ * gaussWeight(k);
            KEE     = KEE     + KEE_GP     * detJ * gaussWeight(k);
            KEtaT   = KEtaT   + KEtaT_GP   * detJ * gaussWeight(k);
            KEtaEta = KEtaEta + KEtaEta_GP * detJ * gaussWeight(k);
            KCX     = KCX     + KCX_GP     * detJ * gaussWeight(k);
            KCC     = KCC     + KCC_GP     * detJ * gaussWeight(k);
            KGX     = KGX     + KGX_GP     * detJ * gaussWeight(k);
            KGC     = KGC     + KGC_GP     * detJ * gaussWeight(k);
            KGG     = KGG     + KGG_GP     * detJ * gaussWeight(k);
            KJX     = KJX     + KJX_GP     * detJ * gaussWeight(k);
            KJG     = KJG     + KJG_GP     * detJ * gaussWeight(k);
            KJJ     = KJJ     + KJJ_GP     * detJ * gaussWeight(k);
            
        % energy and energy difference
            elementEnergy.internalEnergy = elementEnergy.internalEnergy + internalEnergyN1 * detJ * gaussWeight(k);
            elementEnergy.deltaU = elementEnergy.deltaU + (internalEnergyN1 - internalEnergyN) * detJ * gaussWeight(k);
        % entropy and entropy difference
            elementEnergy.entropy = elementEnergy.entropy + entropyN1 * detJ * gaussWeight(k);
            elementEnergy.deltaS = elementEnergy.deltaS + (entropyN1 - entropyN) * detJ * gaussWeight(k);
            


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
    rData{3} = RE;
    rData{4} = REta;
    rData{5} = RC;
    rData{6} = RG;
    rData{7} = RJ;
    kData{1, 1} = KXX;   kData{1, 2} = KXT;   kData{1, 3} = KXE;   kData{1, 4} = KXEta;   kData{1, 5} = KXC;   kData{1, 6} = KXG;   kData{1, 7} = KXJ;   
    kData{2, 1} = KTX;   kData{2, 2} = KTT;   kData{2, 3} = KTE;   kData{2, 4} = KTEta;   kData{2, 5} = KTC;   kData{2, 6} = KTG;   kData{2, 7} = KTJ;
    kData{3, 1} = KEX;   kData{3, 2} = KET;   kData{3, 3} = KEE;   kData{3, 4} = KEEta;   kData{3, 5} = KEC;   kData{3, 6} = KEG;   kData{3, 7} = KEJ;
    kData{4, 1} = KEtaX; kData{4, 2} = KEtaT; kData{4, 3} = KEtaE; kData{4, 4} = KEtaEta; kData{4, 5} = KEtaC; kData{4, 6} = KEtaG; kData{4, 7} = KEtaJ;
    kData{5, 1} = KCX;   kData{5, 2} = KCT;   kData{5, 3} = KCE;   kData{5, 4} = KCEta;   kData{5, 5} = KCC;   kData{5, 6} = KCG;   kData{5, 7} = KCJ;   
    kData{6, 1} = KGX;   kData{6, 2} = KGT;   kData{6, 3} = KGE;   kData{6, 4} = KGEta;   kData{6, 5} = KGC;   kData{6, 6} = KGG;   kData{6, 7} = KGJ;
    kData{7, 1} = KJX;   kData{7, 2} = KJT;   kData{7, 3} = KJE;   kData{7, 4} = KJEta;   kData{7, 5} = KJC;   kData{7, 6} = KJG;   kData{7, 7} = KJJ;
end
0;
end

function Av = matrixToVoigtUnsymmetric(A)
%% Converts nonsymmetric matrix into voigt notation
    Av = [A(1,1); A(2,2); A(3,3); A(1,2); A(2,3); A(1,3); A(2,1); A(3,2); A(3,1)];
end

function A = voigtToMatrixUnsymmetric(Av)
%% Converts nonsymmetric matrix into voigt notation
    A = [Av(1) Av(4) Av(6); 
         Av(7) Av(2) Av(5); 
         Av(9) Av(8) Av(3)];
end

function [C_o_C_sym_voigt] = sym_voigt(A,B, stressStrainA, stressStrainB)
    A11 = A(1,1); A12 = A(1,2); A13 = A(1,3);
    A21 = A(2,1); A22 = A(2,2); A23 = A(2,3);
    A31 = A(3,1); A32 = A(3,2); A33 = A(3,3);

    if strcmpi(stressStrainA,'strain')
        A12 = 2*A12;
        A21 = 2*A21;
        A13 = 2*A13;
        A31 = 2*A31;
        A23 = 2*A23;
        A32 = 2*A32;
    end

    B11 = B(1,1); B12 = B(1,2); B13 = B(1,3);
    B21 = B(2,1); B22 = B(2,2); B23 = B(2,3);
    B31 = B(3,1); B32 = B(3,2); B33 = B(3,3);

    if strcmpi(stressStrainB,'strain')
        B12 = 2*B12;
        B21 = 2*B21;
        B13 = 2*B13;
        B31 = 2*B31;
        B23 = 2*B23;
        B32 = 2*B32;
    end

    C_o_C_sym_voigt = [ A11*B11   A12*B12   A13*B13  0.5*(A11*B12 + A12*B11)     0.5*(A12*B13 + A13*B12)     0.5*(A11*B13 + A13*B11);
                        A21*B21   A22*B22   A23*B23  0.5*(A21*B22 + A22*B21)     0.5*(A22*B23 + A23*B22)     0.5*(A21*B23 + A23*B21);
                        A31*B31   A32*B32   A33*B33  0.5*(A31*B32 + A32*B31)     0.5*(A32*B33 + A33*B32)     0.5*(A31*B33 + A33*B31);
                        A11*B21   A12*B22   A13*B23  0.5*(A11*B22 + A12*B21)     0.5*(A12*B23 + A13*B22)     0.5*(A11*B23 + A13*B21);
                        A21*B31   A22*B32   A23*B33  0.5*(A21*B32 + A22*B31)     0.5*(A22*B33 + A23*B32)     0.5*(A21*B33 + A23*B31);
                        A11*B31   A12*B32   A13*B33  0.5*(A11*B32 + A12*B31)     0.5*(A12*B33 + A13*B32)     0.5*(A11*B33 + A13*B31)];
end

function [DC] = DiffSec(C)
    C11 = C(1,1); C22 = C(2,2); C33 = C(3,3);
    C12 = C(1,2); C13 = C(1,3); C23 = C(2,3);
    
    DC = [   0        C33     C22      0         -C23        0       ;   
             C33      0       C11      0          0         -C13     ;
             C22      C11     0       -C12        0          0       ;
             0        0      -C12     -1/2*C33    1/2*C13    1/2*C23 ;
            -C23      0       0        1/2*C13   -1/2*C11    1/2*C12 ;
             0       -C13     0        1/2*C23    1/2*C12   -1/2*C22 ];
end

function D = secDiffOperator(A)
D = zeros(6, 6);
D(1, 2) = A(3, 3);
D(1, 3) = A(2, 2);
D(1, 5) = -2 * A(2, 3);
D(2, 1) = A(3, 3);
D(2, 3) = A(1, 1);
D(2, 6) = -2 * A(1, 3);
D(3, 1) = A(2, 2);
D(3, 2) = A(1, 1);
D(3, 4) = -2 * A(1, 2);
D(4, 3) = -2 * A(1, 2);
D(4, 4) = -2 * A(3, 3);
D(4, 5) = 2 * A(1, 3);
D(4, 6) = 2 * A(2, 3);
D(5, 1) = -2 * A(2, 3);
D(5, 4) = 2 * A(1, 3);
D(5, 5) = -2 * A(1, 1);
D(5, 6) = 2 * A(1, 2);
D(6, 2) = -2 * A(1, 3);
D(6, 4) = 2 * A(2, 3);
D(6, 5) = 2 * A(1, 2);
D(6, 6) = -2 * A(2, 2);
end

function [C_o_C_special_voigt] = special_voigt(A,B)
    A11 = A(1,1); A12 = A(1,2); A13 = A(1,3);
    A21 = A(2,1); A22 = A(2,2); A23 = A(2,3);
    A31 = A(3,1); A32 = A(3,2); A33 = A(3,3);

    B11 = B(1,1); B12 = B(1,2); B13 = B(1,3);
    B21 = B(2,1); B22 = B(2,2); B23 = B(2,3);
    B31 = B(3,1); B32 = B(3,2); B33 = B(3,3);

    C_o_C_special_voigt = [ A11*B11     A12*B21     A13*B31     0.5*(A11*B21 + A12*B11)   0.5*(A12*B31 + A13*B21)   0.5*(A11*B31 + A13*B11)   0.5*(A11*B21 + A12*B11)   0.5*(A12*B31 + A13*B21)   0.5*(A11*B31 + A13*B11);
                            A21*B12     A22*B22     A23*B32     0.5*(A21*B22 + A22*B12)   0.5*(A22*B32 + A23*B22)   0.5*(A21*B32 + A23*B12)   0.5*(A21*B22 + A22*B12)   0.5*(A22*B32 + A23*B22)   0.5*(A21*B32 + A23*B12);
                            A31*B13     A32*B23     A33*B33     0.5*(A31*B23 + A32*B13)   0.5*(A32*B33 + A33*B23)   0.5*(A31*B33 + A33*B13)   0.5*(A31*B23 + A32*B13)   0.5*(A32*B33 + A33*B23)   0.5*(A31*B33 + A33*B13);
                            A11*B12     A12*B22     A13*B32     0.5*(A11*B22 + A12*B12)   0.5*(A12*B32 + A13*B22)   0.5*(A11*B32 + A13*B12)   0.5*(A11*B22 + A12*B12)   0.5*(A12*B32 + A13*B22)   0.5*(A11*B32 + A13*B12);
                            A21*B13     A22*B23     A23*B33     0.5*(A21*B23 + A22*B13)   0.5*(A22*B33 + A23*B23)   0.5*(A21*B33 + A23*B13)   0.5*(A21*B23 + A22*B13)   0.5*(A22*B33 + A23*B23)   0.5*(A21*B33 + A23*B13);
                            A11*B13     A12*B23     A13*B33     0.5*(A11*B23 + A12*B13)   0.5*(A12*B33 + A13*B23)   0.5*(A11*B33 + A13*B13)   0.5*(A11*B23 + A12*B13)   0.5*(A12*B33 + A13*B23)   0.5*(A11*B33 + A13*B13);
                            A21*B11     A22*B21     A23*B31     0.5*(A21*B21 + A22*B11)   0.5*(A22*B31 + A23*B21)   0.5*(A21*B31 + A23*B11)   0.5*(A21*B21 + A22*B11)   0.5*(A22*B31 + A23*B21)   0.5*(A21*B31 + A23*B11);
                            A31*B12     A32*B22     A33*B32     0.5*(A31*B22 + A32*B12)   0.5*(A32*B32 + A33*B22)   0.5*(A31*B32 + A33*B12)   0.5*(A31*B22 + A32*B12)   0.5*(A32*B32 + A33*B22)   0.5*(A31*B32 + A33*B12);
                            A31*B11     A32*B21     A33*B31     0.5*(A31*B21 + A32*B11)   0.5*(A32*B31 + A33*B21)   0.5*(A31*B31 + A33*B11)   0.5*(A31*B21 + A32*B11)   0.5*(A32*B31 + A33*B21)   0.5*(A31*B31 + A33*B11)];
end

function [A_obar_BVoigt] = A_obar_BVoigt(A,B)
    A_obar_BVoigt = zeros(3,3,3,3);
    map.VoigtUnsymmetric = [1, 5, 9, 4, 8, 7, 2, 6, 3]';
    map.VoigtUnsymmetricTensor = [1,1; 2,2; 3,3; 1,2; 2,3; 1,3; 2,1; 3,2; 3,1];
    Av = A(map.VoigtUnsymmetric);
    Bv = B(map.VoigtUnsymmetric);
    IJ = 0;
    for I = 1:3
        for J = 1:3
            IJ = IJ + 1;
 %           KL = 0;
            for K = 1:3
                for L = 1:3
%                    KL = KL + 1;
                    A_obar_BVoigt(I,J,K,L) = A(I,K)*B(L,J);
                end
            end
        end
    end 
    sub2ind(size(C), [1 2], [1 2], [1 2], [1 2])
end
function C = OuterProduct(A, B)
    C = reshape(A(:) * B(:).', [size(A), size(B)]);
end

function [BQ] = BMatrixTemp(A, b)
    BQ = zeros(6,8);
    Ax = A(1,:);
    Ay = A(2,:);
    Az = A(3,:);
    b1 = b(1);
    b2 = b(2);
    b3 = b(3);
    BQ(1,:) = Ax * b1;
    BQ(2,:) = Ay * b2;
    BQ(3,:) = Az * b3;
    BQ(4,:) = (Ax * b2 + Ay * b1);
    BQ(5,:) = (Ay * b3 + Az * b2);
    BQ(6,:) = (Ax * b3 + Az * b1);
end