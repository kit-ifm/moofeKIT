function Re = gaussDisplacementSCSaintVenantEndpoint(edN1,k,dimension,J,dNrAll,DMat,gaussWeight,I)
indx = dimension*k-(dimension-1):dimension*k;
dNX = (J(:,indx)')\dNrAll(indx,:);
numberOfNodes = size(dNX,2);
edN1 = reshape(edN1,dimension,numberOfNodes);
detJ = det(J(:,indx)');
if detJ < 10*eps
    error('Jacobi determinant equal or less than zero.')
end
FAkt = edN1*dNX';
% B-matrix
numberOfDisplacementDofs = numberOfNodes*dimension;
numberOfSymmetricVoigtDofs = dimension + (dimension^2-dimension)/2;
BAkt = zeros(numberOfSymmetricVoigtDofs,numberOfDisplacementDofs);
BAkt(1,1:dimension:end) = FAkt(1,1)*dNX(1,:);
BAkt(1,2:dimension:end) = FAkt(2,1)*dNX(1,:);
BAkt(1,3:dimension:end) = FAkt(3,1)*dNX(1,:);
BAkt(2,1:dimension:end) = FAkt(1,2)*dNX(2,:);
BAkt(2,2:dimension:end) = FAkt(2,2)*dNX(2,:);
BAkt(2,3:dimension:end) = FAkt(3,2)*dNX(2,:);
BAkt(3,1:dimension:end) = FAkt(1,3)*dNX(3,:);
BAkt(3,2:dimension:end) = FAkt(2,3)*dNX(3,:);
BAkt(3,3:dimension:end) = FAkt(3,3)*dNX(3,:);
BAkt(4,1:dimension:end) = FAkt(1,1)*dNX(2,:) + FAkt(1,2)*dNX(1,:);
BAkt(4,2:dimension:end) = FAkt(2,1)*dNX(2,:) + FAkt(2,2)*dNX(1,:);
BAkt(4,3:dimension:end) = FAkt(3,1)*dNX(2,:) + FAkt(3,2)*dNX(1,:);
BAkt(5,1:dimension:end) = FAkt(1,2)*dNX(3,:) + FAkt(1,3)*dNX(2,:);
BAkt(5,2:dimension:end) = FAkt(2,2)*dNX(3,:) + FAkt(2,3)*dNX(2,:);
BAkt(5,3:dimension:end) = FAkt(3,2)*dNX(3,:) + FAkt(3,3)*dNX(2,:);
BAkt(6,1:dimension:end) = FAkt(1,1)*dNX(3,:) + FAkt(1,3)*dNX(1,:);
BAkt(6,2:dimension:end) = FAkt(2,1)*dNX(3,:) + FAkt(2,3)*dNX(1,:);
BAkt(6,3:dimension:end) = FAkt(3,1)*dNX(3,:) + FAkt(3,3)*dNX(1,:);% Cauchy-Green tensor
CAkt = FAkt'*FAkt;
% Green-Lagrange tensor
En1   = 0.5*(CAkt-I);
En1_v = [En1(1,1) En1(2,2) En1(3,3) 2*En1(1,2) 2*En1(3,2) 2*En1(3,1)]';
% Stresses
DW1_v = 1/2*DMat*En1_v;
% Residual
Re = 2*BAkt'*DW1_v*detJ*gaussWeight(k);
% Tangent
% DW1 = [DW1_v(1) DW1_v(4) DW1_v(6); DW1_v(4) DW1_v(2) DW1_v(5); DW1_v(6) DW1_v(5) DW1_v(3)];
% D2W1 = 1/4*DMat;
% A1 = 2*dNX'*DW1*dNX*detJ*gaussWeight(k);
% MAT = zeros(numberOfDOFs);
% for g = 1:dimension
%     MAT(g:dimension:numberOfDOFs,g:dimension:numberOfDOFs) = A1;
% end
% Ke = 4*BAkt'*D2W1*BAkt*detJ*gaussWeight(k) + MAT;
end