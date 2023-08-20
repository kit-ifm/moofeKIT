function [rData, kData, elementEnergy, array] = deadLoad(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
%DEADLOAD Computes residual resulting from 'dead load' Neumann boundaries

% load objects
shapeFunctionObject = obj.shapeFunctionObject;
meshObject = obj.meshObject;

% loadGeometry
loadGeometry = obj.loadGeometry;
assert(any(strcmp(loadGeometry, {'line', 'area'})), ['Not implemented yet for given loadGeometry (loadGeometry=', loadGeometry, ')!']);

% edof
edof = meshObject.edof(e, :);

% load vector
timeEvaluationPoint = setupObject.time;
if strcmp(setupObject.integrator, {'ExplicitEuler'})
    timeEvaluationPoint = timeEvaluationPoint - setupObject.timeStepSize;
elseif any(strcmp(setupObject.integrator, {'Midpoint', 'DiscreteGradient'}))
    timeEvaluationPoint = timeEvaluationPoint - 1 / 2 * setupObject.timeStepSize;
end
loadVector = obj.loadVector * obj.timeFunction(timeEvaluationPoint);

% gauss points & shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

% nodal coordinates
numberOfDisplacementDofs = size(obj.masterObject.meshObject.nodes, 2) - obj.masterObject.additionalFields;
edR = obj.masterObject.meshObject.nodes(edof, 1:numberOfDisplacementDofs).';
if numberOfDisplacementDofs ~= 3
    edR = [edR; zeros(3-numberOfDisplacementDofs, size(edR, 2))];
end

% projection type
switch obj.projectionType
    case 'none'
    case 'x'
        edR = edR(1, :);
    case 'y'
        edR = edR(2, :);
    case 'z'
        edR = edR(3, :);
    otherwise
        error('type of projection not implemented')
end

% Jacobian
JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);

% initialize residual
RX = rData{1};

% Gauss loop
for k = 1:numberOfGausspoints
    J = JAll(:, size(dN_xi_k_I, 1)*(k - 1)+1:size(dN_xi_k_I, 1)*k);
    if strcmp(loadGeometry, 'line')
        detJ = norm(J);
    elseif strcmp(loadGeometry, 'area')
        detJ = norm(cross(J(:, 1), J(:, 2)));
    end
    RX = RX - kron(N_k_I(k, :).', loadVector) * detJ * gaussWeight(k);
end

% element energy
displacementDofsPerNode = size(obj.masterObject.qR, 2) - obj.masterObject.additionalFields;
fieldsToConsider = 1:displacementDofsPerNode;
if strcmp(obj.loadPhysics, 'thermal')
    if isa(obj.masterObject, 'solidThermoClass')
        fieldsToConsider = displacementDofsPerNode + 1;
    elseif isa(obj.masterObject, 'solidElectroThermoClass')
        fieldsToConsider = displacementDofsPerNode + 2;
    else
        error('Not implemented yet!');
    end
elseif strcmp(obj.loadPhysics, 'electrical')
    if isa(obj.masterObject, 'solidElectroClass') || isa(obj.masterObject, 'solidElectroThermoClass')
        fieldsToConsider = displacementDofsPerNode + 1;
    else
        error('Not implemented yet!');
    end
end
uN1 = (obj.masterObject.qN1(edof, fieldsToConsider) - obj.masterObject.qR(edof, fieldsToConsider)).'; % displacement of nodes
F = -RX; % force on nodes
elementEnergy.externalEnergy = uN1(:).' * F;

% Pass computation data
if ~computePostData
    rData{1} = RX;
end

end
