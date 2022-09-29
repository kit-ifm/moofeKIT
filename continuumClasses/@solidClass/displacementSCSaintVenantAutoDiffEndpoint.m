function displacementSCSaintVenantAutoDiffEndpoint(obj,setupObject)
globalFullEdof = obj.globalFullEdof;
edof = obj.edof;
qN1 = obj.qN1;
numberOfElements = size(globalFullEdof,1);
numberOfDOFs = size(globalFullEdof,2);
dimension = obj.dimension;
flagNumericalTangent = obj.flagNumericalTangent;
strainEnergy = 0;
initializeDataFE(obj);
for e = 1:numberOfElements
    Re = zeros(numberOfDOFs,1);
    RXinitial = zeros(numberOfDOFs,1);
    Ke = zeros(numberOfDOFs);
    %
    mechanicalDOFs = 1:numberOfDOFs;
    numberOfMechanicalDOFs = numel(mechanicalDOFs);
    numberOfFields = 1;     % X
    dofsData = cell(numberOfFields,1);
    dofsData{1,1} = numberOfMechanicalDOFs;
    rDataInitial = cell(numberOfFields,1);
    kDataInitial = cell(numberOfFields,numberOfFields);
    for ii = 1:numberOfFields
        rDataInitial{ii,1} = zeros(dofsData{ii,1},1);
        for jj = 1:numberOfFields
            kDataInitial{ii,jj} = zeros(dofsData{ii,1},dofsData{jj,1});
        end
    end
    %% residual and tangent computation
%     [rData,kData,strainEnergy] = Residuum(e,obj,rDataInitial,kDataInitial,strainEnergy,[],false);

    edof = obj.edof;

    %% Shape functions
    numberOfGausspoints = obj.numberOfGausspoints;
    gaussWeight = obj.shapeFunctions.gaussWeight;
    dimension = obj.dimension;
    dNrAll = obj.shapeFunctions.dNr;

    %
    if flagNumericalTangent
        edN1 = dofsNT.edN1;
    else
        edN1  = obj.qN1(edof(e,:),1:dimension)';
    end
    qR = obj.qR;

    %% Create residual and tangent
    lambda = obj.materialObject.lambda;
    mu = obj.materialObject.mu;
    DMat = [lambda+2*mu lambda lambda 0 0 0;...
        lambda lambda+2*mu lambda 0 0 0;...
        lambda lambda lambda+2*mu 0 0 0;...
        0 0 0 mu 0 0;...
        0 0 0 0 mu 0;...
        0 0 0 0 0 mu];
    I = eye(dimension);


    %% Initialization
    J = qR(edof(e,:),1:dimension)'*dNrAll';
    % Run through all Gauss points
    for k = 1:numberOfGausspoints
%         edN1Dlarray = dlarray(edN1,'SSCB');
        edN1Dlarray = dlarray(edN1(:));
%         [RX,autoDiffGradient] = dlfeval(@ResiduumGaussLevel,edN1Dlarray,obj,dNrAll,J,k,RXinitial,DMat,gaussWeight);
        [RX] = admDiffVFor(@ResiduumGaussLevel,edN1Dlarray,obj,dNrAll,J,k,RXinitial,DMat,gaussWeight);
    end

    Re(mechanicalDOFs,1) = RX;
    Ke(mechanicalDOFs,mechanicalDOFs) = KXX;
    storeDataFE(obj,Re,Ke,globalFullEdof,e);
end
obj.ePot(setupObject.timeStep).strainEnergy = strainEnergy;
end

function [RX,autoDiffGradient] = ResiduumGaussLevel2(edN1,obj,dNrAll,J,k,RX,DMat,gaussWeight)
edN1 = reshape(edN1,3,8);
dimension = obj.dimension;
indx = dimension*k-(dimension-1):dimension*k;
detJ = det(J(:,indx)');
if detJ < 10*eps
    error('Jacobi determinant equal or less than zero.')
end
dNx = (J(:,indx)')\dNrAll(indx,:);
FAkt = edN1*dNx';
% B-matrix
BAkt = BMatrix(dNx,FAkt);
% Cauchy-Green tensor
CAkt = FAkt'*FAkt;
% Green-Lagrange tensor
En1   = 0.5*(CAkt-eye(3));
En1_v = [En1(1,1) En1(2,2) En1(3,3) 2*En1(1,2) 2*En1(3,2) 2*En1(3,1)]';
% Stresses
DW1_v = 1/2*DMat*En1_v;
% DW1 = [DW1_v(1) DW1_v(4) DW1_v(6); DW1_v(4) DW1_v(2) DW1_v(5); DW1_v(6) DW1_v(5) DW1_v(3)];
% Residual
RX = RX + 2*BAkt'*DW1_v*detJ*gaussWeight(k);
% % Tangent
% D2W1 = 1/4*DMat;
% A1 = 2*dNx'*DW1*dNx*detJ*gaussWeight(k);
% MAT = zeros(numberOfDOFs);
% for g = 1:dimension
%     MAT(g:dimension:numberOfDOFs,g:dimension:numberOfDOFs) = A1;
% end
% KXX = KXX + 4*BAkt'*D2W1*BAkt*detJ*gaussWeight(k) + MAT;
% % 
% RX = RX(1);
% % 
autoDiffGradient = dlgradient(RX,edN1(:));
end

function [RX] = ResiduumGaussLevel(edN1,obj,dNrAll,J,k,RX,DMat,gaussWeight)
edN1 = reshape(edN1,3,8);
dimension = obj.dimension;
indx = dimension*k-(dimension-1):dimension*k;
detJ = det(J(:,indx)');
if detJ < 10*eps
    error('Jacobi determinant equal or less than zero.')
end
dNx = (J(:,indx)')\dNrAll(indx,:);
FAkt = edN1*dNx';
% B-matrix
BAkt = BMatrix(dNx,FAkt);
% Cauchy-Green tensor
CAkt = FAkt'*FAkt;
% Green-Lagrange tensor
En1   = 0.5*(CAkt-eye(3));
En1_v = [En1(1,1) En1(2,2) En1(3,3) 2*En1(1,2) 2*En1(3,2) 2*En1(3,1)]';
% Stresses
DW1_v = 1/2*DMat*En1_v;
% DW1 = [DW1_v(1) DW1_v(4) DW1_v(6); DW1_v(4) DW1_v(2) DW1_v(5); DW1_v(6) DW1_v(5) DW1_v(3)];
% Residual
RX = RX + 2*BAkt'*DW1_v*detJ*gaussWeight(k);
end