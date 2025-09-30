function [rData, kData, elementEnergy, array] = mixedSCpHCGJMooneyRivlinModifiedVol2LinearImplicit(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
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

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
mixedFEObject = obj.mixedFEObject;
meshObject = obj.meshObject;
DT = setupObject.timeStepSize;

% aquire general data
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
M_k_I = mixedFEObject.shapeFunctionObject.N_k_I;

numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussWeight = shapeFunctionObject.gaussWeight;

edof = meshObject.edof(e, :);
numberOfXDofs = size(meshObject.globalFullEdof, 2) - size(mixedFEObject.globalEdof, 2);
numberOfInternalNodes = size(M_k_I, 2);

dimension = obj.dimension;

% aquire material data
a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;
% 
alpha = materialObject.alpha;
beta = materialObject.beta;
gamma = materialObject.gamma;
epsilon1 = materialObject.epsilon1;
epsilon2 = materialObject.epsilon2;
% 
rho = materialObject.rhoLinearImplicit;

% aquire the nodal values of the variables for the current element
% phi
edR = obj.qR(edof, 1:dimension).';
edN = obj.qN(edof, 1:dimension).';
phiN = edN(:);
edN1 = dofs.edN1;
phiN1 = edN1(:);
% v
edvN1 = dofs.vN1;
vN1 = edvN1(:);
edvN = obj.qN(edof, dimension+1:end).';
vN = edvN(:);
% midpoint evaluation of the velocity
edvN05 = 0.5*(edvN + edvN1); %(edN1-edN)/DT;
vN05 = edvN05(:);
% update
edHeadN05 = edN + DT/2*edvN; % edN05 = 0.5*(edN+edN1);
phiHeadN05 = edHeadN05(:);

% initialize elementEnergy
elementEnergy.strainEnergy = 0;
elementEnergy.kineticEnergy = 0;

I = eye(dimension);
numberOfDOFs = dimension*size(edof, 2);

% compute Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

edAlphaN = mixedFEObject.qN(e, :).';
edAlphaN1 = dofs.edAlphaN1.';

extractedCNv = edAlphaN(1:6*numberOfInternalNodes);
extractedGNv = edAlphaN(6*numberOfInternalNodes+1:12*numberOfInternalNodes);
extractedJN = edAlphaN(12*numberOfInternalNodes+1:13*numberOfInternalNodes);
extractedCN1v = edAlphaN1(1:6*numberOfInternalNodes);
extractedGN1v = edAlphaN1(6*numberOfInternalNodes+1:12*numberOfInternalNodes);
extractedJN1 = edAlphaN1(12*numberOfInternalNodes+1:13*numberOfInternalNodes);

% initialize residual
RX = rData{1};
RV = rData{2};
RDotCv = rData{3};
RDotGv = rData{4};
RDotJ = rData{5};

% initialize tangent
% KXX = kData{1, 1};
% KCC = kData{2, 2};
% KCG = kData{2, 3};
% KGC = kData{3, 2};
% KJJ = kData{4, 4};

I6 = eye(6);
IV = [1 1 1 0 0 0];
Mphi = zeros(numberOfDOFs);
Mrho = zeros(numberOfDOFs);
MC = zeros(numberOfInternalNodes*6);
MJ = zeros(numberOfInternalNodes);
Phi = zeros(numberOfDOFs,numberOfInternalNodes*6);
Psi = zeros(numberOfDOFs,numberOfInternalNodes*6);
Sigma = zeros(numberOfDOFs,numberOfInternalNodes);
dU_Cv = zeros(numberOfInternalNodes*6,1);
dU_Gv = zeros(numberOfInternalNodes*6,1);
dU_J = zeros(numberOfInternalNodes,1);
ddU_Cv = zeros(numberOfInternalNodes*6);
ddU_Gv = zeros(numberOfInternalNodes*6);
ddU_J = zeros(numberOfInternalNodes);

% Run through all Gauss points
for k = 1:numberOfGausspoints
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_I = computedN_X_I(dN_xi_k_I, J, k);
    % compute the values of the variables at the current Gauss point
    CNv = reshape(extractedCNv, 6, []) * M_k_I(k, :)';
    Mbar = M_k_I(k,:)';
    Mkron = kron(M_k_I(k,:)',I6);
    % CNv = Mkron'*extractedCNv;
    CN = voigtToMatrix(CNv, 'stress');
    CN1v = reshape(extractedCN1v, 6, []) * M_k_I(k, :)';
    CN1 = voigtToMatrix(CN1v, 'stress');
    CN05 = 0.5 * (CN + CN1);

    GNv = reshape(extractedGNv, 6, []) * M_k_I(k, :)';
    GN = voigtToMatrix(GNv, 'stress');
    GN1v = reshape(extractedGN1v, 6, []) * M_k_I(k, :)';
    GN1 = voigtToMatrix(GN1v, 'stress');
    GN05 = 0.5 * (GN + GN1);

    JN = extractedJN.' * M_k_I(k, :)';
    JN1 = extractedJN1.' * M_k_I(k, :)';
    JN05 = 0.5 * (JN + JN1);

    % Deformation gradient
    FxN1 = edN1 * dN_X_I';
    FxN05 = edHeadN05 * dN_X_I';
    % B-matrix (current configuration)
    BN05 = BMatrix(dN_X_I, FxN05);

    % Right Cauchy-Green tensor
    CxN1 = FxN1.' * FxN1;
    CxN05 = FxN05.' * FxN05;
    % Cofactor
    GxN1 = 0.5 * wedge(CxN1, CxN1);
    GxN05 = 0.5 * wedge(CxN05, CxN05);
    % Determinant
    JxN1 = det(FxN1);
    JxN05 = det(FxN05);
    % Derivative of the strain energy function
    dW_CN1 = alpha * trace(CN1) * I;
    dW_GN1 = beta * trace(GN1) * I;
    dW_C = alpha * trace(CN05) * I;
    dW_G = beta * trace(GN05) * I;
    % dW_J = - gamma*JN05^(-1) + 2*epsilon1*epsilon2*(JN05^(2*epsilon2-1) - JN05^(-2*epsilon2-1));
    % dW_JN1 = - gamma*JN1^(-1) + 2*epsilon1*epsilon2*(JN1^(2*epsilon2-1) - JN1^(-2*epsilon2-1));
    dW_J = c*(JN05 - 1) - d;
    dW_JN1 = c*(JN1 - 1) - d;
    % Second Piola Kirchhoff stress tensor
    SN1 = 2 * (dW_CN1 + wedge(dW_GN1, CxN1) + 1/2*dW_JN1*JxN1^(-1)*GxN1);
    SN05 = 2 * (dW_C + wedge(dW_G, CxN05) + 1/2*dW_J*JxN05^(-1)*GxN05);
    FvN05 = edvN05 * dN_X_I';
    CvN05 = FxN05.'*FvN05 + FvN05.'*FxN05;
    GvN05 = wedge(CxN05,CvN05);
    JvN05 = 1/2*JxN05^(-1)*GxN05(:).'*CvN05(:);
    SN05_v = matrixToVoigt(SN05, 'stress');
    Nkron = kron(N_k_I(k,:)',I);
    if ~computePostData
        % RX = RX + (Nkron*Nkron')*(1/DT*(phiN1-phiHeadN05)-0.5*vN1) * detJ * gaussWeight(k);
        % RV = RV + ((rho*Nkron*Nkron')*1/DT*(vN1-vN) + BN05.' * SN05_v) * detJ * gaussWeight(k);
        % RDotCv = RDotCv + kron(M_k_I(k, :)',matrixToVoigt(1/DT*(CN1 - CN) - CvN05,'stress')) * detJ * gaussWeight(k);
        % RDotGv = RDotGv + kron(M_k_I(k, :)',matrixToVoigt(1/DT*(GN1 - GN) - GvN05,'stress')) * detJ * gaussWeight(k);
        % RDotJ = RDotJ + M_k_I(k, :)' * (1/DT*(JN1 - JN) - JvN05) * detJ * gaussWeight(k);
        % Residual
        Mphi = Mphi + Nkron*Nkron'*detJ*gaussWeight(k);
        Mrho = Mrho + rho*(Nkron*Nkron')*detJ*gaussWeight(k);
        MC = MC + Mkron*Mkron'*detJ*gaussWeight(k);
        MJ = MJ + M_k_I(k, :)'*M_k_I(k, :)*detJ*gaussWeight(k);
        Phi = Phi + 2*BN05'*Mkron'*detJ*gaussWeight(k);
        Psi = Psi + 2*BN05'*secDiffOperator(CxN05)*Mkron'*detJ*gaussWeight(k);
        Sigma = Sigma + BN05.'*JxN05^(-1)*matrixToVoigt(GxN05,'stress')*Mbar.'*detJ*gaussWeight(k);
        dU_Cv = dU_Cv + Mkron*matrixToVoigt(dW_C,'stress')*detJ*gaussWeight(k);
        dU_Gv = dU_Gv + Mkron*matrixToVoigt(dW_G,'stress')*detJ*gaussWeight(k);
        dU_J = dU_J + M_k_I(k, :)'*dW_J*detJ*gaussWeight(k);
        ddU_Cv = ddU_Cv + Mkron*alpha*0.5*kron(Mbar,IV'*IV)'*detJ*gaussWeight(k);
        ddU_Gv = ddU_Gv + Mkron*beta*0.5*kron(Mbar,IV'*IV)'*detJ*gaussWeight(k);
        ddU_J = ddU_J + M_k_I(k, :)'*c*0.5*M_k_I(k, :)*detJ*gaussWeight(k);

        % Strain energy
        W = alpha/2 * ((trace(CN1))^2 - 9) + beta/2 * ((trace(GN1))^2 - 9) + c / 2 * (JN1 - 1)^2 - d * (JN1-1);
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + W * detJ * gaussWeight(k);
        
    else
        % stress at gausspoint
        PN1 = FxN1 * SN1;
        stressTensor.FirstPK = PN1;
        stressTensor.Cauchy = 1 / det(FxN1) * PN1 * FxN1';
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end
if ~computePostData
    % pass residual
    % rData{1} = 1/DT*(phiN1-phiHeadN05)-0.5*vN1;
    % rData{1} = Mphi*(1/DT*(phiN1-phiN)-0.5*(vN + vN1));
    % rData{1} = 1/DT*(phiN1-phiHeadN05)-0.5*vN1;
    rData{1} = Mphi*(1/DT*(phiN1-phiHeadN05)-0.5*vN1);
    rData{2} = Mrho*1/DT*(vN1-vN) + Phi*(MC\dU_Cv) + Psi*(MC\dU_Gv) + Sigma*(MJ\dU_J);
    rData{3} = MC*1/DT*(extractedCN1v-extractedCNv) - Phi.'*vN05;
    rData{4} = MC*1/DT*(extractedGN1v-extractedGNv) - Psi.'*vN05;
    rData{5} = MJ*1/DT*(extractedJN1-extractedJN) - Sigma.'*vN05;   
    % rData{1} = RX;
    % rData{2} = RV;
    % rData{3} = RDotCv;
    % rData{4} = RDotGv;
    % rData{5} = RDotJ;
    
    % pass tangent
    kData{1, 1} = 1/DT*Mphi;        kData{1, 2} = -0.5*Mphi;
    kData{2, 2} = Mrho*1/DT;        kData{2, 3} = Phi*(MC\ddU_Cv);     kData{2, 4} = Psi*(MC\ddU_Gv);   kData{2, 5} = Sigma*(MJ\ddU_J);
    kData{3, 2} = -Phi.'*0.5;       kData{3, 3} = MC*1/DT;
    kData{4, 2} = -Psi.'*0.5;       kData{4, 4} = MC*1/DT;
    kData{5, 2} = -Sigma.'*0.5;     kData{5, 5} = MJ*1/DT;

    elementEnergy.kineticEnergy = elementEnergy.kineticEnergy + 0.5*vN1'*Mrho*vN1;
    % elementEnergy.kineticEnergy = elementEnergy.kineticEnergy + 0.5*vN1'*Mrho*vN1 + 0.5*extractedCN1v'*MC*extractedCN1v + 0.5*extractedGN1v'*MC*extractedGN1v + 0.5*extractedJN1'*MJ*extractedJN1;
end
end
% function D = secDiffOperator(A)
% D = zeros(6, 6);
% D(1, 2) = A(3, 3);
% D(1, 3) = A(2, 2);
% D(1, 5) = -2*A(2, 3);
% D(2, 1) = A(3, 3);
% D(2, 3) = A(1, 1);
% D(2, 6) = -2*A(1, 3);
% D(3, 1) = A(2, 2);
% D(3, 2) = A(1, 1);
% D(3, 4) = -2*A(1, 2);
% D(4, 3) = -2*A(1, 2);
% D(4, 4) = -2*A(3, 3);
% D(4, 5) = 2*A(1, 3);
% D(4, 6) = 2*A(2, 3);
% D(5, 1) = -2*A(2, 3);
% D(5, 4) = 2*A(1, 3);
% D(5, 5) = -2*A(1, 1);
% D(5, 6) = 2*A(1, 2);
% D(6, 2) = -2*A(1, 3);
% D(6, 4) = 2*A(2, 3);
% D(6, 5) = 2*A(1, 2);
% D(6, 6) = -2*A(2, 2);
% end


% function D = secDiffOperator(A)
% D = zeros(6, 6);
% D(1, 2) = A(3, 3);
% D(1, 3) = A(2, 2);
% D(1, 5) = -0.5 * (A(2, 3) + A(3, 2));
% D(2, 1) = A(3, 3);
% D(2, 3) = A(1, 1);
% D(2, 6) = -0.5 * (A(1, 3) + A(3, 1));
% D(3, 1) = A(2, 2);
% D(3, 2) = A(1, 1);
% D(3, 4) = -0.5 * (A(1, 2) + A(2, 1));
% D(4, 3) = -0.5 * (A(1, 2) + A(2, 1));
% D(4, 4) = -0.5 * A(3, 3);
% D(4, 5) = 0.25 * (A(1, 3) + A(3, 1));
% D(4, 6) = 0.25 * (A(2, 3) + A(3, 2));
% D(5, 1) = -0.5 * (A(2, 3) + A(3, 2));
% D(5, 4) = 0.25 * (A(1, 3) + A(3, 1));
% D(5, 5) = -0.5 * A(1, 1);
% D(5, 6) = 0.25 * (A(1, 2) + A(2, 1));
% D(6, 2) = -0.5 * (A(1, 3) + A(3, 1));
% D(6, 4) = 0.25 * (A(2, 3) + A(3, 2));
% D(6, 5) = 0.25 * (A(1, 2) + A(2, 1));
% D(6, 6) = -0.5 * A(2, 2);
% end

function D = secDiffOperator(A)
D = zeros(6, 6);
D(1, 2) = A(3, 3);
D(1, 3) = A(2, 2);
D(1, 5) = -2 * (A(2, 3) + A(3, 2));
D(2, 1) = A(3, 3);
D(2, 3) = A(1, 1);
D(2, 6) = -2 * (A(1, 3) + A(3, 1));
D(3, 1) = A(2, 2);
D(3, 2) = A(1, 1);
D(3, 4) = -2 * (A(1, 2) + A(2, 1));
D(4, 3) = -(A(1, 2) + A(2, 1));
D(4, 4) = -A(3, 3);
D(4, 5) = (A(1, 3) + A(3, 1));
D(4, 6) = (A(2, 3) + A(3, 2));
D(5, 1) = -(A(2, 3) + A(3, 2));
D(5, 4) = (A(1, 3) + A(3, 1));
D(5, 5) = -A(1, 1);
D(5, 6) = (A(1, 2) + A(2, 1));
D(6, 2) = -(A(1, 3) + A(3, 1));
D(6, 4) = (A(2, 3) + A(3, 2));
D(6, 5) = (A(1, 2) + A(2, 1));
D(6, 6) = -A(2, 2);
end