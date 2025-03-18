%% Pure Bending plate test
% This test allows to test the fullfillment of the Kirchhoff condition
% (gamma = 0) for arbitrary mesh geometries. The Kirchhoff condition is
% satisfied, if the transverse forces (q1, q2) vanish for the implemented
% case of pure bending.
% Fullfilled by: ANS elements
% Not fullfilled by: EAS element

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'plateSimple';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.postPlotType = 'stress';
setupObject.plotObject.stress = struct('type', 'Cauchy', 'component', 13);
setupObject.newton.tolerance = 1e-7;

dofObject = dofClass; % required object for dof and object handling

%% continuum Objects
meshType = 'regular'; % regular, distortedLine, distortedPatch
plateLength = 50;
s = 5;
orderShapeFunctions = 1;
plateObject = plateClass(dofObject);

if strcmp(meshType, 'regular')
    [plateObject.meshObject.nodes, plateObject.meshObject.edof, boundaryEdof] = meshRectangle(plateLength, plateLength, 2, 2, orderShapeFunctions);
    plateObject.meshObject.nodes = plateObject.meshObject.nodes + plateLength / 2;
elseif strcmp(meshType, 'distortedLine')
    [plateObject.meshObject.nodes, plateObject.meshObject.edof, boundaryEdof] = meshRectangle(plateLength, plateLength, 2, 1, orderShapeFunctions);
    nodeToDistort1 = find(plateObject.meshObject.nodes(:, 1) == 0 & plateObject.meshObject.nodes(:, 2) == -plateLength/2);
    plateObject.meshObject.nodes(nodeToDistort1, 1) = plateObject.meshObject.nodes(nodeToDistort1, 1) + s;
    nodeToDistort2 = find(plateObject.meshObject.nodes(:, 1) == 0 & plateObject.meshObject.nodes(:, 2) == plateLength/2);
    plateObject.meshObject.nodes(nodeToDistort2, 1) = plateObject.meshObject.nodes(nodeToDistort2, 1) - s;
    plateObject.meshObject.nodes = plateObject.meshObject.nodes + plateLength / 2;
elseif strcmp(meshType, 'distortedPatch')
    [plateObject.meshObject.nodes, plateObject.meshObject.edof, boundaryEdof] = meshPureBendingDistorted(plateLength, 4, orderShapeFunctions, false);
end


plateObject.materialObject.name = 'Hooke';
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 1e5;
plateObject.materialObject.nu = 0;
plateObject.h = 1;
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = orderShapeFunctions;
plateObject.shapeFunctionObject.numberOfGausspoints = (orderShapeFunctions + 1)^2;
plateObject.elementDisplacementType = 'displacement'; % displacement, eas, selectiveReducedIntegration
plateObject.elementNameAdditionalSpecification = 'PetrovGalerkinBatheDvorkin'; % PetrovGalerkinBatheDvorkin, BatheDvorkin, SimoRifai
plateObject.mixedFEObject.condensation = true;
plateObject.mixedFEObject.typeShapeFunctionData = 4;

% Dirichlet boundary
boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(plateObject.meshObject.nodes(:, 1) == 0);
boundary1.nodalDof = [1, 3];
boundary1.masterObject = plateObject;

boundary2 = dirichletClass(dofObject);
boundary2.nodeList = find(plateObject.meshObject.nodes(:, 1) == plateLength);
boundary2.nodalDof = [1, 3];
boundary2.masterObject = plateObject;

boundary3 = dirichletClass(dofObject);
boundary3.nodeList = find(plateObject.meshObject.nodes(:, 1) == 0);
boundary3.nodalDof = 2;
boundary3.timeFunction = @(t, z) 1;
boundary3.masterObject = plateObject;

boundary4 = dirichletClass(dofObject);
boundary4.nodeList = find(plateObject.meshObject.nodes(:, 1) == plateLength);
boundary4.nodalDof = 2;
boundary4.timeFunction = @(t, z) -1;
boundary4.masterObject = plateObject;

%% solver
dofObject = runNewton(setupObject, dofObject);
disp(['Maximum vertical displacement: ', num2str(max(abs(plateObject.qN1(:, 1))))]);

function [nodes, edof, bounEdof] = meshPureBendingDistorted(length, numberOfElementsPerDirection, orderShapeFunctions, serendipity)
bounEdof = [];
edof = [];
nodes = [];

for ii = 1:numberOfElementsPerDirection
    for jj = 1:numberOfElementsPerDirection
        [nodesElem, edofElem] = meshPatchTestDistorted2D(length/numberOfElementsPerDirection, length/numberOfElementsPerDirection, orderShapeFunctions, serendipity);
        nodesElem(:, 1) = nodesElem(:, 1) + length / numberOfElementsPerDirection * (ii - 1);
        nodesElem(:, 2) = nodesElem(:, 2) + length / numberOfElementsPerDirection * (jj - 1);
        if ~isempty(edof)
            edofElem = edofElem + max(max(edof));
            [nodes, edof] = meshTie(nodes, edof, nodesElem, edofElem, eps);
        else
            nodes = nodesElem;
            edof = edofElem;
        end
    end
end
end