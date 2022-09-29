function dataFE = displacementHookeEndpoint(obj,setupObject)
%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = linearEndPoint(obj,'PropertyName',PropertyValue)
%
% Description
%
% homogenous linear-elastic isotropic strain-energy function, evaluated at time n+1, i.e. implicid euler method.
%
% 08.02.2015 M.FRANKE

%% Check input

%% Shape functions
globalFullEdof = obj.globalFullEdof;
edof = obj.edof;
numberOfGausspoints = obj.numberOfGausspoints;
gaussWeight = obj.shapeFunctions.gaussWeight;
NAll = obj.shapeFunctions.N';
dNrAll = obj.shapeFunctions.dNr';
qR = obj.qR;
qN1 = obj.qN1;
numberOfElements = size(globalFullEdof,1);
numberOfDOFs = size(globalFullEdof,2);
dimension = obj.dimension;

lambda = obj.materialData.lambda;
mu = obj.materialData.mu;
DMat = [lambda+2*mu lambda lambda 0 0 0;...
    lambda lambda+2*mu lambda 0 0 0;...
    lambda lambda lambda+2*mu 0 0 0;...
    0 0 0 mu 0 0;...
    0 0 0 0 mu 0;...
    0 0 0 0 0 mu];

%% Create residual and tangent
dataFE(numberOfElements) = struct('edofE',[],'Re',[],'ePot',[],'pI',[],'pJ',[],'pK',[]);
parfor e = 1:numberOfElements
    [edofH1,edofH2] = expandEdof(globalFullEdof(e,:));
    dataFE(e).pI = edofH1';
    dataFE(e).pJ = edofH2';
    dataFE(e).edofE = double(globalFullEdof(e,:))';
    % Element routine
    Re = zeros(numberOfDOFs,1);
    Ke = zeros(numberOfDOFs);
    ePot = 0;
    edN1 = qN1(edof(e,:),1:dimension)';
    edRef = qR(edof(e,:),1:dimension)';
    uN1 = edN1 - edRef;
    J = qR(edof(e,:),1:dimension)'*dNrAll;
    % Run through all Gauss points
    for k = 1:numberOfGausspoints
        indx = dimension*k-(dimension-1):dimension*k;
        detJ = det(J(:,indx)');
        if detJ < 10*eps
            error('Jacobi determinant equal or less than zero.')
        end
        % FIXME
%         dNx = (J(:,indx)')\dNrAll(indx,:);
        dNx = (J(:,indx)')\(dNrAll(:,indx))';
        B = zeros(6,numberOfDOFs);
        B(1,1:3:end) = dNx(1,:);
        B(2,2:3:end) = dNx(2,:);
        B(3,3:3:end) = dNx(3,:);
        B(4,1:3:end) = dNx(2,:);
        B(4,2:3:end) = dNx(1,:);
        B(5,2:3:end) = dNx(3,:);
        B(5,3:3:end) = dNx(2,:);
        B(6,1:3:end) = dNx(3,:);
        B(6,3:3:end) = dNx(1,:);
        Re = Re + (B'*DMat*B)*uN1(:)*detJ*gaussWeight(k);
        Ke = Ke + B'*DMat*B*detJ*gaussWeight(k);
    end
    dataFE(e).Re = Re;
    dataFE(e).pK = Ke(:);
    dataFE(e).ePot = ePot;
end
end