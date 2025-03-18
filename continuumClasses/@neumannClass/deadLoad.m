function [rData, kData, elementEnergy, array] = deadLoad(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
%DEADLOAD Computes residual resulting from 'dead load' Neumann boundaries

% load objects
shapeFunctionObject = obj.shapeFunctionObject;
meshObject = obj.meshObject;
masterObject = obj.masterObject;

% edof
edof = meshObject.edof(e, :);
numberOfNodes = length(edof);

% nodal coordinates
numberOfAdditionalDofs = sum(masterObject.dofsPerAdditionalField);
numberOfDisplacementDofs = size(masterObject.qR, 2) - numberOfAdditionalDofs;
numberOfSpatialCoordinates = size(masterObject.meshObject.nodes, 2) - numberOfAdditionalDofs;
edR = masterObject.meshObject.nodes(edof, 1:numberOfSpatialCoordinates).';
if numberOfSpatialCoordinates ~= 3
    edR = [edR; zeros(3-numberOfSpatialCoordinates, numberOfNodes)];
end

% loadGeometry
loadGeometry = obj.loadGeometry;
assert(any(strcmp(loadGeometry, {'line', 'area'})), ['Not implemented yet for given loadGeometry (loadGeometry=', loadGeometry, ')!']);

% loadPhysics
loadPhysics = obj.loadPhysics;
fieldsToConsider = 1:numberOfDisplacementDofs;
if strcmp(loadPhysics, 'thermal')
    if isa(masterObject, 'solidThermoClass')
        fieldsToConsider = numberOfDisplacementDofs + 1;
    elseif isa(masterObject, 'solidElectroThermoClass')
        fieldsToConsider = numberOfDisplacementDofs + 2;
    else
        error('Not implemented yet!');
    end
elseif strcmp(loadPhysics, 'electrical')
    if isa(masterObject, 'solidElectroClass') || isa(masterObject, 'solidElectroThermoClass')
        fieldsToConsider = numberOfDisplacementDofs + 1;
    else
        error('Not implemented yet!');
    end
end

% load vector
timeEvaluationPoint = setupObject.time;
if strcmp(setupObject.integrator, {'ExplicitEuler'})
    timeEvaluationPoint = timeEvaluationPoint - setupObject.timeStepSize;
elseif any(strcmp(setupObject.integrator, {'Midpoint', 'DiscreteGradient'}))
    timeEvaluationPoint = timeEvaluationPoint - 1 / 2 * setupObject.timeStepSize;
end
if ~isempty(obj.loadVector)
    % when loadVector is given
    loadVectorMatrix = kron(eye(numberOfNodes), obj.loadVector * obj.timeFunction(timeEvaluationPoint));
else
    % when loadVectorFunction is given
    loadVectorMatrix = zeros(length(fieldsToConsider)*numberOfNodes, numberOfNodes);
    for i=1:numberOfNodes
        indx = length(fieldsToConsider)*(i-1)+1:length(fieldsToConsider)*(i-1)+length(fieldsToConsider);
        loadVectorMatrix(indx, i) = obj.loadVectorFunction(edR(:, i)) * obj.timeFunction(timeEvaluationPoint);
    end
end

% gauss points & shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
N_k_I = shapeFunctionObject.N_k_I;
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;

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

% check, if object is axisymmetric
isAxisymmetric = false;
if isa(obj.masterObject, 'axisymmetricSolidClass')
    isAxisymmetric = true;
end

% Gauss loop
for k = 1:numberOfGausspoints
    J = JAll(:, size(dN_xi_k_I, 1)*(k - 1)+1:size(dN_xi_k_I, 1)*k);
    if strcmp(loadGeometry, 'line')
        detJ = norm(J);
    elseif strcmp(loadGeometry, 'area')
        detJ = norm(cross(J(:, 1), J(:, 2)));
    end
    if isAxisymmetric
        r = N_k_I(k, :) * edR(1, :).';
        RX = RX - 2 * pi * loadVectorMatrix * N_k_I(k, :).' * r * detJ * gaussWeight(k);
    else
        RX = RX - loadVectorMatrix * N_k_I(k, :).' * detJ * gaussWeight(k);
    end
end

% element energy
uN1 = (obj.masterObject.qN1(edof, fieldsToConsider) - obj.masterObject.qR(edof, fieldsToConsider)).'; % displacement of nodes
F = -RX; % force on nodes
elementEnergy.externalEnergy = uN1(:).' * F;

% Pass computation data
if ~computePostData
    rData{1} = RX;
end

end
