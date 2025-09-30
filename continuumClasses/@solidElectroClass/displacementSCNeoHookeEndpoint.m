function displacementSCNeoHookeEndpoint(obj,setupObject)
%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = disp_neoHooke_endPoint(obj,'PropertyName',PropertyValue)
%
% Description
%
%
% 06.02.2017 M.Franke A.Janz

%% Check input

globalFullEdof = obj.globalFullEdof;
edof = obj.edof;
numberOfGausspoints = obj.numberOfGausspoints;
gaussWeight = obj.shapeFunctions.gaussWeight;
% NAll = obj.shapeFunctions.N;
dNrAll = obj.shapeFunctions.dNr;
qR = obj.qR;
% qN = obj.qN;
qN1 = obj.qN1;
numberOfElements = size(globalFullEdof,1);
numberOfDOFs = size(globalFullEdof,2);
numberOfMechanicalDOFs = numberOfDOFs-8; % mechanical dofs
dimension = obj.dimension;

%% Create residual and tangent
lambda = obj.materialObject.lambda;
mu = obj.materialObject.mu;
c1 = obj.materialObject.c1;
c2 = obj.materialObject.c2;
epsilon0 = obj.materialObject.epsilon0;
I = eye(dimension);
voigt=[1 5 9 4 8 7 2 6 3]'; % Reformulation to voigt notation
shape=[1 7 9 4 2 8 6 5 3]'; % Reformulation to matrix notation (useful for reshape operation)
% permutation matrix
index = [1 34 67 101 134 167 201 234 267 301 334 367 401 434 467 501 534 567 601 634 667 701 734 767 772 808 844 880 916 952 988 1024];
P = sparse(32,32);
P(index) = 1;
flagNumericalTangent = false;
% element loop
omegaEnergy = 0;
storageFEObject = obj.storageFEObject;
initializeDataFE(storageFEObject);
for e = 1:numberOfElements
    % Element routine
    edN1 = qN1(edof(e,:),1:dimension)';
    phiN1 = qN1(edof(e,:),dimension+1)';
    J = qR(edof(e,:),1:dimension)'*dNrAll';
    % gauss loop
    [Re, Ke, omegaEnergy] =  gauss(omegaEnergy,numberOfDOFs,numberOfMechanicalDOFs,numberOfGausspoints,dimension,J,dNrAll,edN1,phiN1,mu,lambda,c1,c2,epsilon0,I,voigt,shape,gaussWeight);
    % numerical tangent
    if flagNumericalTangent
        [Ke] = numericalTangent(omegaEnergy,numberOfDOFs,numberOfMechanicalDOFs,numberOfGausspoints,dimension,J,dNrAll,edN1,phiN1,mu,lambda,c1,c2,epsilon0,I,voigt,shape,gaussWeight,Re);
    end
    Re = P*Re;
    Ke = P*Ke*P';
    storeDataFE(storageFEObject,Re,Ke(:),globalFullEdof,e);
end
obj.elementData(setupObject.timeStep).omegaEnergy = omegaEnergy;
end

function [Re, Ke, omegaEnergy] = gauss(omegaEnergy,numDOFs,numberOfMechanicalDOFs,NGP,DIM,J,dNr,edN1,phiN1,mu,lambda,c1,c2,epsilon0,I,voigt,shape,wp)
Re = zeros(numDOFs,1);
Ke = zeros(numDOFs);
for k = 1:NGP
    indx = DIM*k-(DIM-1):DIM*k;
    detJ = det(J(:,indx)');
    if detJ < 10*eps
        error('Jacobi determinant equal or less than zero.')
    end
    dNx = (J(:,indx)')\dNr(indx,:);
    FN1 = edN1*dNx';
    %
    EN1 = -(phiN1*dNx')';
    % B-matrix displacement
    BN1 = zeros(6,numberOfMechanicalDOFs);
    % deltaC11
    BN1(1,1:3:end) = 2*FN1(1,1)*dNx(1,:);
    BN1(1,2:3:end) = 2*FN1(2,1)*dNx(1,:);
    BN1(1,3:3:end) = 2*FN1(3,1)*dNx(1,:);
    % deltaC22
    BN1(2,1:3:end) = 2*FN1(1,2)*dNx(2,:);
    BN1(2,2:3:end) = 2*FN1(2,2)*dNx(2,:);
    BN1(2,3:3:end) = 2*FN1(3,2)*dNx(2,:);
    % deltaC33
    BN1(3,1:3:end) = 2*FN1(1,3)*dNx(3,:);
    BN1(3,2:3:end) = 2*FN1(2,3)*dNx(3,:);
    BN1(3,3:3:end) = 2*FN1(3,3)*dNx(3,:);
    % 2*deltaC12
    BN1(4,1:3:end) = 2*(FN1(1,1)*dNx(2,:)+FN1(1,2)*dNx(1,:));
    BN1(4,2:3:end) = 2*(FN1(2,1)*dNx(2,:)+FN1(2,2)*dNx(1,:));
    BN1(4,3:3:end) = 2*(FN1(3,1)*dNx(2,:)+FN1(3,2)*dNx(1,:));
    % 2*deltaC23
    BN1(5,1:3:end) = 2*(FN1(1,2)*dNx(3,:)+FN1(1,3)*dNx(2,:));
    BN1(5,2:3:end) = 2*(FN1(2,2)*dNx(3,:)+FN1(2,3)*dNx(2,:));
    BN1(5,3:3:end) = 2*(FN1(3,2)*dNx(3,:)+FN1(3,3)*dNx(2,:));
    % 2*deltaC13
    BN1(6,1:3:end) = 2*(FN1(1,1)*dNx(3,:)+FN1(1,3)*dNx(1,:));
    BN1(6,2:3:end) = 2*(FN1(2,1)*dNx(3,:)+FN1(2,3)*dNx(1,:));
    BN1(6,3:3:end) = 2*(FN1(3,1)*dNx(3,:)+FN1(3,3)*dNx(1,:));
    numberOfElectricalDOFs = numDOFs - numberOfMechanicalDOFs;
    % B-matrix electrical field
    BEN1 = zeros(6,numberOfElectricalDOFs);
    % E1*deltaE1
    BEN1(1,1:end) = -2*EN1(1)*dNx(1,:);
    % E2*deltaE2
    BEN1(2,1:end) = -2*EN1(2)*dNx(2,:);
    % E3*deltaE3
    BEN1(3,1:end) = -2*EN1(3)*dNx(3,:);
    % E2*deltaE1 + E1*deltaE2
    BEN1(4,1:end) = -(EN1(2)*dNx(1,:) + EN1(1)*dNx(2,:));
    % E3*deltaE2 + E2*deltaE3
    BEN1(5,1:end) = -(EN1(3)*dNx(2,:) + EN1(2)*dNx(3,:));
    % E3*deltaE1 + E1*deltaE3
    BEN1(6,1:end) = -(EN1(3)*dNx(1,:) + EN1(1)*dNx(3,:));
    % Cauchy-Green tensor
    CN1 = FN1'*FN1;
    % Inverse strain tensor
    CInvN1 = CN1\I;
    % Invariant
    detFN1 = (det(CN1))^0.5;
    % Strain energy function
    EN1OEN1 = EN1*EN1';
    Omega = mu/2*(trace(CN1)-3) + lambda/2*(log(detFN1))^2 - mu*log(detFN1) + c1*EN1'*EN1 + c2*CN1(voigt)'*EN1OEN1(voigt);
    omegaEnergy = omegaEnergy + Omega*detJ*wp(k);
    % First derivative of strain energy density function
    DCOmegaN1 = mu/2*(I-CInvN1) + lambda/2*log(detFN1)*CInvN1 + c2*EN1OEN1;
    % Voigt notation
    DCOmegaN1_v = [ DCOmegaN1(1,1);...
        DCOmegaN1(2,2);...
        DCOmegaN1(3,3);...
        DCOmegaN1(1,2);...
        DCOmegaN1(2,3);...
        DCOmegaN1(1,3)];
    DEOmegaN1 = 2*c1*EN1 + 2*c2*CN1*EN1;
    % Residual
    Re = Re + [ BN1'*DCOmegaN1_v;...
        -dNx'*DEOmegaN1]*detJ*wp(k);
    % Second derivative of strain energy function
    CInvN1Sym = [CInvN1.*CInvN1 CInvN1.*CInvN1(:,[2 3 1]); CInvN1.*CInvN1([2 3 1],:) (CInvN1.*CInvN1([2 3 1],[2 3 1])+CInvN1(:,[2 3 1]).*CInvN1([2 3 1],:))/2];
    CInvN1_v = [    CInvN1(1,1);...
        CInvN1(2,2);...
        CInvN1(3,3);...
        CInvN1(1,2);...
        CInvN1(2,3);...
        CInvN1(1,3)];
    D2W1 = (mu/2-lambda/2*log(detFN1))*CInvN1Sym+lambda/4*(CInvN1_v*CInvN1_v');
    % Tangent
    A1 = 2*dNx'*DCOmegaN1*dNx*detJ*wp(k);
    MAT = zeros(numberOfMechanicalDOFs);
    for g = 1:DIM
        MAT(g:DIM:numberOfMechanicalDOFs,g:DIM:numberOfMechanicalDOFs) = A1;
    end
    Ke(1:numberOfMechanicalDOFs,1:numberOfMechanicalDOFs) = Ke(1:numberOfMechanicalDOFs,1:numberOfMechanicalDOFs) + BN1'*D2W1*BN1*detJ*wp(k) + MAT;
    Ke(1:numberOfMechanicalDOFs,numberOfMechanicalDOFs+1:end) = Ke(1:numberOfMechanicalDOFs,numberOfMechanicalDOFs+1:end) + c2*BN1'*BEN1*detJ*wp(k);
    Ke(numberOfMechanicalDOFs+1:end,1:numberOfMechanicalDOFs) = Ke(numberOfMechanicalDOFs+1:end,1:numberOfMechanicalDOFs) + c2*BEN1'*BN1*detJ*wp(k);
    Ke(numberOfMechanicalDOFs+1:end,numberOfMechanicalDOFs+1:end) = Ke(numberOfMechanicalDOFs+1:end,numberOfMechanicalDOFs+1:end) + (2*c1*dNx'*dNx + 2*c2*dNx'*CN1*dNx)*detJ*wp(k);
end
end

function [Ke] = numericalTangent(numDOFs,numberOfMechanicalDOFs,NGP,DIM,J,dNr,edN1,phiN1,mu,lambda,c1,c2,epsilon0,I,voigt,shape,wp,Re);
h = 1e-5;
Ke = zeros(numDOFs);
for jj = 1:numDOFs
    edN1Num = edN1;
    phiN1Num = phiN1;
    if jj <= 24
        edN1Num(jj) = edN1Num(jj) + h;
    else
        phiN1Num(jj-24) = phiN1Num(jj-24) + h;
    end
    [ReNum, ~, ~] =  gauss(numDOFs,numberOfMechanicalDOFs,NGP,DIM,J,dNr,edN1Num,phiN1Num,mu,lambda,c1,c2,epsilon0,I,voigt,shape,wp);
    Ke(:,jj) =  Ke(:,jj) + (ReNum - Re)/h;
end
end