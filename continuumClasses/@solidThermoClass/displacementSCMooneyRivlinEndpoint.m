function [rData, kData, elementEnergy, array] = displacementSCMooneyRivlinEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% displacementSCMooneyRivlinEndpoint(obj,setupObject)
%
% Description
% dispCascadeSC=: displacement based cascPade formulation
% based on a symmetric formulation in 2nd PK S and right Cauchy-Green C.
% Mooney Rivlin coupled mechanical and thermal strain-energy function,
% evaluated at time n+1 (implicit euler scheme).
% 05.04.2017 M.S., A.J., M.F.


%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dimension = obj.dimension;
edof = meshObject.edof;
DT = setupObject.timeStepSize;

a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;
kappa = materialObject.kappa;
beta = materialObject.beta;
thetaR = materialObject.thetaR;
k0 = materialObject.k0;

% Initialize sub-residual vector
RX = rData{1};
RT = rData{2};
% Initialize sub-tangent matrices
KXX = kData{1, 1};
KXT = kData{1, 2};
KTX = kData{2, 1};
KTT = kData{2, 2};
% degrees of freedom
qR = obj.qR;
qN = obj.qN;
edR = qR(edof(e,:),1:dimension).';
edN = qN(edof(e,:),1:dimension).';
edN1 = dofs.edN1;
thetaN1 = dofs.thetaN1;
thetaN = qN(edof(e,:),dimension+1).';

elementEnergy.helmholtzEnergy = 0;

JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

    % Temperature
    thetaN1e = N_k_I(k, :) * thetaN1.';
    thetaNe = N_k_I(k, :) * thetaN.';
    
    % Deformation Gradient
    FxN1 = edN1*dN_X_I';
    FxN = edN*dN_X_I';
    CxN1 = FxN1.'*FxN1;
    CxN = FxN.'*FxN;
    if ~numel(CxN1)==1
        GxN1 = 0.5*wedge(CxN1,CxN1);
    else
        GxN1 = 0.5*CxN1*CxN1;
    end
    cxN1 = det(CxN1);    
    cxN = det(CxN);    
    
    % B-matrix (current configuration)
    BN1 = BMatrix(dN_X_I,FxN1);
    BN1 = 2*BN1;

    % Helmholtz energy function
    PsiIsoN1 = a*(trace(CxN1)-dimension) + b*(trace(GxN1)-dimension);
    PsiVolN1 = -d*log(sqrt(cxN1)) + c/2*(sqrt(cxN1)-1)^2;
    PsiThermoN1 = kappa*(thetaN1e - thetaR - thetaN1e*log(thetaN1e/thetaR));
    PsiCoupledN1 = -dimension*beta*(thetaN1e - thetaR)*(c*(sqrt(cxN1) - 1) - d/sqrt(cxN1));
    PsiN1 = PsiIsoN1 + PsiVolN1 + PsiThermoN1 + PsiCoupledN1;
    elementEnergy.helmholtzEnergy = elementEnergy.helmholtzEnergy + PsiN1*detJ*gaussWeight(k);

    %% Impuls balance
    % Derivative of the strain energy function
    DPsi_C = a*eye(dimension);
    DPsi_G = b*eye(dimension);
    DPsi_cxN1 = -d/(2*cxN1) + c/2*(1-1/sqrt(cxN1)) - dimension*beta*(thetaN1e - thetaR)*(c/2*1/sqrt(cxN1) + d/2*(cxN1)^(-3/2));

    if ~numel(CxN1)==1
        SN1 = 2*(DPsi_C + wedge(DPsi_G,CxN1) + DPsi_cxN1*GxN1);
    else
        SN1 = 2*(DPsi_C + DPsi_G*CxN1 + DPsi_cxN1*GxN1);
    end
    
    if ~computePostData
        if ~(numel(CxN1)==1)
            rX = BN1.'*0.5*matrixToVoigt(SN1, 'stress');
        else
            rX = BN1.'*0.5*SN1;
        end
        etaN1 = kappa*log(thetaN1e/thetaR)+dimension*beta*(c*(sqrt(cxN1)-1)-d/sqrt(cxN1));
        etaN = kappa*log(thetaNe/thetaR)+dimension*beta*(c*(sqrt(cxN)-1)-d/sqrt(cxN));
        Q = -k0*(1/cxN1*GxN1*(dN_X_I*thetaN1'));
        rT = (N_k_I(k, :)'*thetaN1e*(etaN1-etaN)/DT-dN_X_I'*Q);
    
        %% summation of tangent and residual
        RX = RX + rX*detJ*gaussWeight(k);
        RT = RT + rT*detJ*gaussWeight(k);

        %% TODO: compute tangent
    else
        %% Stress computation
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
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
