function out = displacementSaintVenantEndpoint(obj,configObject)
%% Creates the residual and the tangent of the given obj.
%
% Syntax
%
% out = stVenantEndPoint(obj,'PropertyName',PropertyValue)
%
% Description
%
% St. Venant Kirchhoff strain-energy function, evaluated at time n+1, i.e. implicid euler method.
%
% 10.01.2012 C.HESCH

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

%% Create residual and tangent
lambda = obj.materialData.lambda;
mu = obj.materialData.mu;
DMat = [lambda+2*mu lambda lambda 0 0 0;...
    lambda lambda+2*mu lambda 0 0 0;...
    lambda lambda lambda+2*mu 0 0 0;...
    0 0 0 mu 0 0;...
    0 0 0 0 mu 0;...
    0 0 0 0 0 mu];
I = eye(dimension);
out(numberOfElements) = struct('edofE',[],'Re',[],'ePot',[],'pI',[],'pJ',[],'pK',[]);
parfor e = 1:numberOfElements
    [edofH1,edofH2] = expandEdof(globalFullEdof(e,:));
    out(e).pI = edofH1';
    out(e).pJ = edofH2';
    out(e).edofE = double(globalFullEdof(e,:))';
    % Element routine
    Re = zeros(numberOfDOFs,1);
    Ke = zeros(numberOfDOFs);
    ePot = 0;
    edN1 = qN1(edof(e,:),1:dimension)';
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
        FAkt = edN1*dNx';
        % B-matrix
        BAkt = zeros(6,numberOfDOFs);
        BAkt(1,1:3:end)=FAkt(1,1)*dNx(1,:);
        BAkt(1,2:3:end)=FAkt(2,1)*dNx(1,:);
        BAkt(1,3:3:end)=FAkt(3,1)*dNx(1,:);
        BAkt(2,1:3:end)=FAkt(1,2)*dNx(2,:);
        BAkt(2,2:3:end)=FAkt(2,2)*dNx(2,:);
        BAkt(2,3:3:end)=FAkt(3,2)*dNx(2,:);
        BAkt(3,1:3:end)=FAkt(1,3)*dNx(3,:);
        BAkt(3,2:3:end)=FAkt(2,3)*dNx(3,:);
        BAkt(3,3:3:end)=FAkt(3,3)*dNx(3,:);
        BAkt(4,1:3:end)=FAkt(1,1)*dNx(2,:)+FAkt(1,2)*dNx(1,:);
        BAkt(4,2:3:end)=FAkt(2,1)*dNx(2,:)+FAkt(2,2)*dNx(1,:);
        BAkt(4,3:3:end)=FAkt(3,1)*dNx(2,:)+FAkt(3,2)*dNx(1,:);
        BAkt(5,1:3:end)=FAkt(1,2)*dNx(3,:)+FAkt(1,3)*dNx(2,:);
        BAkt(5,2:3:end)=FAkt(2,2)*dNx(3,:)+FAkt(2,3)*dNx(2,:);
        BAkt(5,3:3:end)=FAkt(3,2)*dNx(3,:)+FAkt(3,3)*dNx(2,:);
        BAkt(6,1:3:end)=FAkt(1,1)*dNx(3,:)+FAkt(1,3)*dNx(1,:);
        BAkt(6,2:3:end)=FAkt(2,1)*dNx(3,:)+FAkt(2,3)*dNx(1,:);
        BAkt(6,3:3:end)=FAkt(3,1)*dNx(3,:)+FAkt(3,3)*dNx(1,:);
        % Cauchy-Green tensor
        CAkt = FAkt'*FAkt;
        % Green-Lagrange tensor
        En1   = 0.5*(CAkt-I);
        En1_v = [En1(1,1) En1(2,2) En1(3,3) 2*En1(1,2) 2*En1(3,2) 2*En1(3,1)]';
        % Strain energy
        ePot = ePot + 0.5*En1_v'*DMat*En1_v*detJ*gaussWeight(k);
        % Stresses
        DW1_v = 1/2*DMat*En1_v;
        DW1 = [DW1_v(1) DW1_v(4) DW1_v(6); DW1_v(4) DW1_v(2) DW1_v(5); DW1_v(6) DW1_v(5) DW1_v(3)];
        % Residual
        Re = Re + 2*BAkt'*DW1_v*detJ*gaussWeight(k);
        % Tangent
        D2W1 = 1/4*DMat;
        A1 = 2*dNx'*DW1*dNx*detJ*gaussWeight(k);
        MAT = zeros(numberOfDOFs);
        for g = 1:dimension
            MAT(g:dimension:numberOfDOFs,g:dimension:numberOfDOFs) = A1;
        end
        Ke = Ke + 4*BAkt'*D2W1*BAkt*detJ*gaussWeight(k) + MAT;
    end
    out(e).Re = Re;
    out(e).pK = Ke(:);
    out(e).ePot = ePot;
end
end