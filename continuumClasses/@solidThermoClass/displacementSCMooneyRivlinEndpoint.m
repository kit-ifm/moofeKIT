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


map.Voigt = [1 5 9 4 8 7]';
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;
N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
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
edN1 = dofs.edN1;
edN = qN(edof(e,:),1:dimension).';
thetaN1 = dofs.thetaN1;
thetaN = qN(edof(e,:),dimension+1).';

elementEnergy.helmholtzEnergy = 0;

J = qR(edof(e,:),1:dimension)'*dNr';
JN1 = edN1 * dNr';
% Run through all Gauss points
for k = 1:numberOfGausspoints
    index = dimension*k-(dimension-1):dimension*k;
    detJ = det(J(:,index)');
    detJN1 = det(JN1(:, index)');
    if detJ < 10*eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:, index)') \ dNr(index, :);
    % Temperature
    thetaN1e = N(k, :) * thetaN1.';
    thetaNe = N(k, :) * thetaN.';
    
    % Deformation Gradient
    FxN1 = edN1*dNx';
    FxN = edN*dNx';
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
    BN1 = BMatrix(dNx,FxN1);
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
        Q = -k0*(1/cxN1*GxN1*(dNx*thetaN1'));
        rT = (N(k, :)'*thetaN1e*(etaN1-etaN)/DT-dNx'*Q);
    
        %% summation of tangent and residual
        RX = RX + rX*detJ*gaussWeight(k);
        RT = RT + rT*detJ*gaussWeight(k);

        %% TODO: compute tangent
    else
        %% Stress computation
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N, k, gaussWeight, detJ, detJN1, stressTensor, setupObject, dimension);
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
