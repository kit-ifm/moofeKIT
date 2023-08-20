function [rData, kData, elementEnergy, array] = displacementNeoHookeViscoMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)

%% Creates the residual and the tangent of the given obj.
%
% Description:
%   material model: viscoelastic extention of Neo Hookian  material model
%       -> lambda and mu describe hyperelastic part
%       -> lambda_Visco and mu_Visco desribe viscoelastic part
%   viscoelasticity via internal variable C_i (inelastic part of Cauchy
%     Green strain tensor)
%
% 25.01.2023 Moritz Hille: first version of viscoelastic material behavior

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;

% element degree of freedom tables and more
edof = meshObject.edof(e, :);
dimension = obj.dimension;
DT = setupObject.timeStepSize;
dimensionVisco = sum(1:dimension); % number of DOF internal variable on GP -> dimension of local residual

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

%initialize energy
% TODO

%% extract dofs for element
% nodal dofs
qR = obj.qR;
qN = obj.qN;
qN1 = obj.qN1;
edR = qR(edof, 1:dimension)';
edN = qN(edof, 1:dimension)';
edN1 = qN1(edof, 1:dimension)';

%% material data
lambda = materialObject.lambda;
mu = materialObject.mu;
lambdaVisco = materialObject.lambdaVisco; % viscous Lamé-parameter lambda_e
muVisco = materialObject.muVisco; % viscous Lamé-parameter mu_e
vVol = 500;
vDev = 100;
V1 = 1 / (2 * vDev);
V2 = 1 / (dimension^2 * vVol) - 1 / (2 * dimension * vDev);

%% initialize residual and tangent (global)
elementEnergy.internalEnergy = 0;
elementEnergy.deltaU = 0;

RG = rData{1}; % G = global = Residual for external (displacement) dofs
KGG = kData{1, 1};
dimR = numel(rData{1});

% assignment displacement dofs
correlation = zeros(dimR/dimension, 2); % auxiliary matrix for assignment of displacement dof for numerical tangent
for mm = 1:dimR / dimension
    correlation((3 * mm)-2, 1) = mm;
    correlation((3 * mm)-1, 1) = mm;
    correlation(3*mm, 1) = mm;
    correlation((3 * mm)-2, 2) = 1;
    correlation((3 * mm)-1, 2) = 2;
    correlation(3*mm, 2) = 3;
end

%% Run through all Gauss points
for k = 1:numberOfGausspoints
    % initialize residual and tangent (local) for additional internal variables
    RL = zeros(dimension, dimension); % L = local  = Residual for internal (viscous) dofs
    KLL = zeros(dimensionVisco);
    KLG = zeros(dimensionVisco, dimR);
    KGL = zeros(dimR, dimensionVisco);
    KGG_GP = zeros(dimR); % needed for static condensation

    % Jacobian matrix and determinant
    [detJ, detJStruct, dN_X_I, ~, ~, ~] = computeAllJacobian(edR, edN, edN1, dN, k, setupObject);

    % deformation gradient
    edN05 = 1 / 2 * (edN + edN1);
    FxN05 = edN05 * dN_X_I';
    CxN05 = FxN05.' * FxN05;


    % internal element dofs
    %     obj.epsilonViscoN1(e,k,:) = [1, 1, 1, 0, 0, 0];
    CiN1V(:, :) = obj.epsilonViscoN(e, k, :);

    if ~computePostData

        numTangEps = 1i * 1e-50; % increment for complex numerical tangent
        CiNV = CiN1V;

        %% energy
        %         elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (B * uN1(:))' * C * (B * uN1(:)) * detJ * gaussWeight(k);

        %% internal Newton's method
        setupObject.newtonVisco.step = 0;
        while setupObject.newtonVisco.step <= setupObject.newtonVisco.maximumSteps
            setupObject.newtonVisco.step = setupObject.newtonVisco.step + 1;

            % case 0 = computation residual, case 1 to dimensionVisco = computation column ii of numerical tangent
            for ii = 0:dimensionVisco
                CiN1Vakt = CiN1V;
                % inelastic Cauchy Green
                if ii > 0
                    CiN1Vakt(ii) = CiN1Vakt(ii) + numTangEps;
                end
                CiN1 = voigtToMatrix(CiN1Vakt, 'strain');
                CiN = voigtToMatrix(CiNV, 'strain');
                CiN05 = 1 / 2 * (CiN + CiN1);
                CiN05inv = CiN05 \ eye(dimension);

                % auxiliary variables
                LambdaN05 = CxN05 / CiN05;
                JViscoN05 = sqrt(det(LambdaN05));

                % first derivatives of internal energy: (e = psiCom + psiVis)
                DPsiVisDLambdaN05 = muVisco / 2 * (eye(3) - eye(3) / (LambdaN05.')) + lambdaVisco / 2 * (log(JViscoN05) + JViscoN05^2 - JViscoN05) * (eye(3) / (LambdaN05.'));
                DPsiVisDCiN05 = -0.5 * ((LambdaN05.') * DPsiVisDLambdaN05 * CiN05inv + CiN05inv * (DPsiVisDLambdaN05.') * LambdaN05);
                if ii == 0
                    DPsiVisDCiN05Resi = DPsiVisDCiN05;
                end
                % stresses
                GammaN05 = -DPsiVisDCiN05; % inelastic stress tensor
                MN05 = 2 * CiN05 * GammaN05; % Mandel stress tensor

                % 4th order viscous resilience tensor (V^-1):MN05
                VViscoN05 = V1 * MN05.' + V2 * trace(MN05) * eye(dimension);

                % local residual and tangent
                if ii == 0 % case ii = 0: compute local residual RL
                    RL = (CiN1 - CiN) / DT - 2 * CiN05 * VViscoN05;
                    RLV = matrixToVoigt(RL, 'strain');
                else % case ii > 0: compute column of local tangent (numerically)
                    RLimag = (CiN1 - CiN) / DT - 2 * CiN05 * VViscoN05;
                    KLL(:, ii) = imag(matrixToVoigt(RLimag, 'strain')) / imag(numTangEps);
                end
            end

            %% update local variable
            RLnorm = RL;
            %             RLnorm = DT*sum(sum(DPsiVisDCiN05Resi.*RL));
            normresi = norm(RLnorm);
            %             fprintf('local residual = %9.6g \n',normresi);
            if normresi < 1e-6
                setupObject.newtonVisco.step = 100;
            elseif setupObject.newtonVisco.step == 40 || normresi > 10^(10)
                error('Newton iteration failed')
            else
                CiN1V = CiN1V - KLL \ RLV;
            end
        end
        obj.epsilonViscoN1(e, k, :) = CiN1V;

        %% global residual, energy and tangent
        % case 0 = computation global residual and internal energy, case 1 to dimR = computation column ii of numerical tangent of KGG
        % case dimR+1 to dimR + dimensionVisco = computation column ii of  numerical tangent of KGL
        for ii = 0:dimR + dimensionVisco
            edN1Akt = edN1;
            CiN1Vakt = CiN1V;
            if (ii > 0 && ii <= dimR)
                jj = correlation(ii, 2);
                kk = correlation(ii, 1);
                edN1Akt(jj, kk) = edN1Akt(jj, kk) + numTangEps;
            elseif ii > dimR
                CiN1Vakt(ii-dimR) = CiN1Vakt(ii-dimR) + numTangEps;
            end

            %modified Cauchy Green and internal viscous variable
            edN05 = 1 / 2 * (edN1Akt + edN);
            FxN05 = edN05 * dN_X_I';
            FxN1 = edN1Akt * dN_X_I';
            FxN = edN * dN_X_I';
            CxN05 = FxN05.' * FxN05;
            CxN1 = FxN1.' * FxN1;
            CxN = FxN.' * FxN;
            CxN05inv = CxN05 \ eye(dimension);
            CiN1 = voigtToMatrix(CiN1Vakt, 'strain');
            CiN = voigtToMatrix(CiNV, 'strain');
            CiN05 = 1 / 2 * (CiN + CiN1);
            CiN05inv = CiN05 \ eye(dimension);

            % B-matrix (midpoint)
            BN05 = BMatrix(dN_X_I, FxN05);

            % auxiliary variables
            JN05 = sqrt(det(CxN05));
            JN1 = sqrt(det(CxN1));
            JN = sqrt(det(CxN));
            LambdaN05 = CxN05 / CiN05;
            LambdaN1 = CxN1 / CiN1;
            LambdaN = CxN / CiN;
            JViscoN05 = sqrt(det(LambdaN05));
            JViscoN1 = sqrt(det(LambdaN1));
            JViscoN = sqrt(det(LambdaN));

            % first derivatives of internal energy: (e(C,Ci) = psiCom(C) + psiVis(C,Ci))
            DPsiComDCxN05 = mu / 2 * (eye(3) - CxN05inv) + lambda / 2 * (log(JN05) + JN05^2 - JN05) * CxN05inv;
            DPsiVisDLambdaN05 = muVisco / 2 * (eye(3) - eye(3) / (LambdaN05.')) + lambdaVisco / 2 * (log(JViscoN05) + JViscoN05^2 - JViscoN05) * (eye(3) / (LambdaN05.'));
            DPsiVisDCxN05 = 0.5 * (DPsiVisDLambdaN05 * CiN05inv + CiN05inv * DPsiVisDLambdaN05.');
            DPsiVisDCiN05 = -0.5 * ((LambdaN05.') * DPsiVisDLambdaN05 * CiN05inv + CiN05inv * (DPsiVisDLambdaN05.') * LambdaN05);


            % stresses
            GammaN05 = -DPsiVisDCiN05; % inelastic stress tensor
            MN05 = 2 * CiN05 * GammaN05; % Mandel stress tensor
            SN05 = 2 * (DPsiComDCxN05 + DPsiVisDCxN05); % 2nd PK stress tensor
            %             SN05 = mu * (eye(3) - CxN05inv) + lambda * CxN05inv * (log(JN05) + JN05^2 - JN05) ...
            %                 +muVisco * (CiN05inv - CxN05inv) + lambdaVisco * CxN05inv * (log(JViscoN05) + JViscoN05^2 - JViscoN05); % 2nd PK stress tensor
            SN05V = matrixToVoigt(SN05, 'stress'); % 2nd PK stress tensor in Voigt notation


            % 4th order viscous resilience tensor (V^-1):MN05
            VViscoN05 = V1 * MN05.' + V2 * trace(MN05) * eye(dimension);

            % global residual, internal energy and remaining tangent parts KGG, KGL & KLG
            if ii == 0 % case ii = 0: compute local residual RL
                RG = RG + BN05.' * SN05V * detJ * gaussWeight(k);
                % internal energy e(C,Ci) = psiCom(C) + psiVis(C,Ci)
                psiVolComN1 = internalEnergyVolumetricPart(JN1, lambda);
                psiVolComN = internalEnergyVolumetricPart(JN, lambda);
                psiVolVisN1 = internalEnergyVolumetricPart(JViscoN1, lambdaVisco);
                psiVolVisN = internalEnergyVolumetricPart(JViscoN, lambdaVisco);
                psiComN1 = mu / 2 * (trace(CxN1) - dimension - 2 * log(JN1)) + psiVolComN1;
                psiComN = mu / 2 * (trace(CxN) - dimension - 2 * log(JN)) + psiVolComN;
                psiVisN1 = muVisco / 2 * (trace(LambdaN1) - dimension - 2 * log(JViscoN1)) + psiVolVisN1;
                psiVisN = muVisco / 2 * (trace(LambdaN) - dimension - 2 * log(JViscoN)) + psiVolVisN;
                elementEnergy.internalEnergy = elementEnergy.internalEnergy + (psiComN1 + psiVisN1) * detJ * gaussWeight(k);
                elementEnergy.deltaU = elementEnergy.deltaU + (psiComN1 - psiComN + psiVisN1 - psiVisN) * gaussWeight(k);
            elseif (ii > 0 && ii <= dimR) % case 0 < ii <= dimR: compute column of KGG & KGL (numerically)
                RGimag = BN05.' * SN05V;
                RLimag = (CiN1 - CiN) / DT - 2 * CiN05 * VViscoN05;
                KGG_GP(:, ii) = imag(RGimag) / imag(numTangEps);
                KLG(:, ii) = imag(matrixToVoigt(RLimag, 'strain')) / imag(numTangEps);
            else % case ii > dimR: compute column of KLG
                RGimag = BN05.' * SN05V;
                RLimag = (CiN1 - CiN) / DT - 2 * CiN05 * VViscoN05;
                KGL(:, ii-dimR) = imag(RGimag) / imag(numTangEps);
                KLL(:, ii-dimR) = imag(matrixToVoigt(RLimag, 'strain')) / imag(numTangEps);
            end
        end

        %% global tangent (static condensation)
        KGG = KGG + (KGG_GP - KGL * (KLL \ KLG)) * detJ * gaussWeight(k);


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
    rData{1} = RG;
    kData{1, 1} = KGG;
end
end

function [psiVol] = internalEnergyVolumetricPart(J, lambda)
psiVol = lambda / 2 * (log(J)^2 + (J - 1)^2);
end
