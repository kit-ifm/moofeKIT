function Re = gaussDisplacementHookeEnpoint(uN1,k,dimension,J,dNrAll,DMat,gaussWeight)
indx = dimension*k-(dimension-1):dimension*k;
detJ = det(J(:,indx)');
if detJ < 10*eps
    error('Jacobi determinant equal or less than zero.')
end
dNx = (J(:,indx)')\dNrAll(indx,:);
numberOfNodes = size(dNx,2);
numberOfDisplacementDofs = numberOfNodes*dimension;
numberOfSymmetricVoigtDofs = dimension + (dimension^2-dimension)/2;
B = zeros(numberOfSymmetricVoigtDofs,numberOfDisplacementDofs);
B(1,1:dimension:end) = dNx(1,:);
B(2,2:dimension:end) = dNx(2,:);
B(3,3:dimension:end) = dNx(3,:);
B(4,1:dimension:end) = dNx(2,:);
B(4,2:dimension:end) = dNx(1,:);
B(5,2:dimension:end) = dNx(3,:);
B(5,3:dimension:end) = dNx(2,:);
B(6,1:dimension:end) = dNx(3,:);
B(6,3:dimension:end) = dNx(1,:);
Re = (B'*DMat*B)*uN1*detJ*gaussWeight(k);
% Ke = (B'*DMat*B)*detJ*gaussWeight(k);
end