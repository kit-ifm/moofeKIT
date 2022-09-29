function [detJ, detJStruct, dN_X_I, J, JN, JN1] = computeAllJacobian(edR,edN,edN1,dN_xi_k_I,k,setupObject)
% COMPUTEJACOBIAN function to compute the Jacobian determinant, derivative
% of the shapefunctions wrt x and the Jacobian matrix
%
% CALL [detJ, detJStruct, index] = computeJacobian(J,JN,JN1,dimension,k,setupObject)
% input:    edR dentoes the the geometric nodes of the considered element for reference configuration
%           edN dentoes the the geometric nodes of the considered element for configuration at tN
%           edN1 dentoes the the geometric nodes of the considered element for configuration at tN1
%           dN_xi_k_I array of 1st derivative of the shapefunctions wrt xi
%           k denotes the k-th Gauss point
%           setupObject is needed for the tolerance of the Jacobian determinant
% output:   detJ Jacobian determinant for reference configuration
%           detJStruct structure with   'R' for Jacobian determinant for reference configuration
%                                       'N' for Jacobian determinant for configuration tN
%                                       'N1' for Jacobian determinant for current configuration tN1
%           dN_X_I array of 1st derivative of the shapefunctions wrt X at gausspoint
%           J, JN, JN1: Jacobi matrix due to reference (R) configuration, configuration tN (N) and current configuration (N1)

tolerance = setupObject.toleranceDetJ;
% dN_xi_I = squeezeMoofeKIT(dN_xi_k_I(:,k,:));
dN_xi_I = reshape(dN_xi_k_I(:,k,:),[size(dN_xi_k_I,1),size(dN_xi_k_I,3)]);
% compute Jacobi matrix and determinant
[J,detJ] = computeJacobian(edR,dN_xi_I,tolerance,setupObject.computePostData);
[JN,detJN] = computeJacobian(edN,dN_xi_I,tolerance,setupObject.computePostData);
[JN1,detJN1] = computeJacobian(edN1,dN_xi_I,tolerance,setupObject.computePostData);
detJStruct = struct('R',detJ,'N',detJN,'N1',detJN1);
% derivative of the shape functions with respect to the physical coodrinate X (or x) at gauss point k
dN_X_I = J' \ dN_xi_I;
end
                 