function array = beamMassMatrixElement(obj, setupObject, e)

%% compute mass matrix for element e
% load objects
meshObject = obj.meshObject;
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
% element degree of freedom tables
edof = meshObject.edof;
globalFullEdof = meshObject.globalFullEdof;
% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
% nodal dofs
qR = meshObject.nodes;
displacement_dofs_per_node = obj.numberOfDofsPerNode;
additionalFields = obj.additionalFields;
% material data and voigt notation
if isfield(materialObject,'rho') && isfield(materialObject,'A')
    rho = materialObject.rho;
    area = materialObject.A;
    rhoA = rho * area;
elseif isfield(materialObject,'rhoA')
    rhoA = materialObject.rhoA;
else
    error("No valid input of inertial parameters")
end

if isfield(materialObject,'rho') && isfield(materialObject,'I')
    rho = materialObject.rho;
    inertia = materialObject.I;
    rhoI = rho * inertia;
elseif isfield(materialObject,'rhoI')
    rhoI = materialObject.rhoI;
else
    error("No valid input of inertial parameters")
end

numberOfDofs = numel(globalFullEdof(e, :));
Me = zeros(numberOfDofs);
edR = qR(edof(e, :), 1:displacement_dofs_per_node)';

if ( strcmpi(obj.theory,'Timoshenko') || strcmpi(obj.theory,'GeometricallyExact'))
    has_rotational_inertia = true;
else
    has_rotational_inertia = false;
end


for k = 1:numberOfGausspoints
    dN_xi_I = reshape(dN_xi_k_I(:, k, :), [size(dN_xi_k_I, 1), size(dN_xi_k_I, 3)]);
    % compute Jacobi matrix and determinant
    [~, detJ] = computeJacobian(edR, dN_xi_I, setupObject.toleranceDetJ);

    A1 = (N_k_I(k, :)' * N_k_I(k, :)) * rhoA;
    MAT = zeros(numberOfDofs);
    sizeN = size(N_k_I, 2) * (displacement_dofs_per_node + additionalFields);
    for l = 1:displacement_dofs_per_node
        MAT(l:displacement_dofs_per_node+additionalFields:sizeN, l:displacement_dofs_per_node+additionalFields:sizeN) = A1;
    end
    
    if has_rotational_inertia
        % additional DOF (rotations) induce inertial effects as well
        A2 = (N_k_I(k, :)' * N_k_I(k, :)) * rhoI;
        MAT((displacement_dofs_per_node+1):displacement_dofs_per_node+additionalFields:sizeN, (displacement_dofs_per_node+1):displacement_dofs_per_node+additionalFields:sizeN) = A2;
    end

    if isa(obj, 'beamVelocityClass')
        % TO CHECK: do sth else here?
    end

    Me = Me + MAT * detJ * gaussWeight(k);

end

array.Me = Me;

end

