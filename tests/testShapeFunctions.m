% unit test for shape functions

disp('=======================================');
disp('Now testing: shape functions');
disp('=======================================');

%%%%%%%%%%
% Settings
%%%%%%%%%%
% General
tolerance = 1e-13;

numberOfGausspoints1D = 2;
numberOfGausspoints2D = 4;
numberOfGausspoints3D = 4;

% Lagrange shape functions
numberOfNodes1D = [1, 2, 3];
numberOfNodes2D = [1, 3, 4, 6, 8, 9];
numberOfNodes3D = [1, 4, 8, 10, 20, 27];

% Metric shape functions
numberOfNodesMetricShapeFunctions1D = [2, 3];
numberOfNodesMetricShapeFunctions2D = [4, 8, 9];
numberOfNodesMetricShapeFunctions3D = [8, 20];

%%%%%%%%%%%%%%%%%%%%%%%%%
% Check and process input
%%%%%%%%%%%%%%%%%%%%%%%%%
assert(size(numberOfNodes1D, 1) == 1 || size(numberOfNodes1D, 2) == 1, 'numberOfNodes1D is dimensioned incorrectly!');
assert(size(numberOfNodes2D, 1) == 1 || size(numberOfNodes2D, 2) == 1, 'numberOfNodes2D is dimensioned incorrectly!');
assert(size(numberOfNodes3D, 1) == 1 || size(numberOfNodes3D, 2) == 1, 'numberOfNodes3D is dimensioned incorrectly!');
assert(size(numberOfGausspoints1D, 1) == 1 || size(numberOfGausspoints1D, 2) == 1, 'numberOfGausspoints1D is dimensioned incorrectly!');
assert(size(numberOfGausspoints2D, 1) == 1 || size(numberOfGausspoints2D, 2) == 1, 'numberOfGausspoints2D is dimensioned incorrectly!');
assert(size(numberOfGausspoints3D, 1) == 1 || size(numberOfGausspoints3D, 2) == 1, 'numberOfGausspoints3D is dimensioned incorrectly!');
assert(size(numberOfNodesMetricShapeFunctions1D, 1) == 1 || size(numberOfNodesMetricShapeFunctions1D, 2) == 1, 'numberOfNodesMetricShapeFunctions1D is dimensioned incorrectly!');
assert(size(numberOfNodesMetricShapeFunctions2D, 1) == 1 || size(numberOfNodesMetricShapeFunctions2D, 2) == 1, 'numberOfNodesMetricShapeFunctions2D is dimensioned incorrectly!');
assert(size(numberOfNodesMetricShapeFunctions3D, 1) == 1 || size(numberOfNodesMetricShapeFunctions3D, 2) == 1, 'numberOfNodesMetricShapeFunctions3D is dimensioned incorrectly!');
assert(all(ismember(numberOfNodesMetricShapeFunctions1D, numberOfNodes1D)), 'numberOfNodesMetricShapeFunctions1D must be a subset of numberOfNodes1D!');
assert(all(ismember(numberOfNodesMetricShapeFunctions2D, numberOfNodes2D)), 'numberOfNodesMetricShapeFunctions2D must be a subset of numberOfNodes2D!');
assert(all(ismember(numberOfNodesMetricShapeFunctions3D, numberOfNodes3D)), 'numberOfNodesMetricShapeFunctions3D must be a subset of numberOfNodes3D!');

numberOfNodesCell = {numberOfNodes1D, numberOfNodes2D, numberOfNodes3D};
numberOfGausspointsCell = {numberOfGausspoints1D, numberOfGausspoints2D, numberOfGausspoints3D};

numberOfNodesMetricShapeFunctionsCell = {numberOfNodesMetricShapeFunctions1D, numberOfNodesMetricShapeFunctions2D, numberOfNodesMetricShapeFunctions3D};

maxNumberOfNodes = max([length(numberOfNodes1D), length(numberOfNodes2D), length(numberOfNodes3D)]);
maxNumberOfGausspoints = max([length(numberOfGausspoints1D), length(numberOfGausspoints2D), length(numberOfGausspoints3D)]);
lagrangeShapeFunctionCell = cell(3, maxNumberOfNodes, maxNumberOfGausspoints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symbolic computation of Lagrange shape functions and their derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for dimension = 1:3
    numberOfNodesArray = numberOfNodesCell{dimension};
    numberOfGausspointsArray = numberOfGausspointsCell{dimension};
    for ii = 1:length(numberOfNodesArray)
        numberOfNodes = numberOfNodesArray(ii);
        elementGeometryType = determineElementGeometryType(dimension, numberOfNodes);
        for jj = 1:length(numberOfGausspointsArray)
            numberOfGausspoints = numberOfGausspointsArray(jj);
            gaussPoints = gaussPointsAndWeights(dimension, numberOfGausspoints, elementGeometryType);
            [N_k_I, dN_xi_k_I, d2N_xi_xi_k_I] = computeLagrangeShapeFunctionSymbolic(dimension, numberOfNodes, numberOfGausspoints, gaussPoints);
            lagrangeShapeFunctionCell{dimension, ii, jj} = struct('N_k_I', N_k_I, 'dN_xi_k_I', dN_xi_k_I, 'd2N_xi_xi_k_I', d2N_xi_xi_k_I);

            % sum of all shape functions must equal 1
            assert(all(sum(N_k_I, 2)-ones(numberOfGausspoints, 1) < tolerance), ['Error in Lagrange shape function (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ', numberOfGausspoints=', num2str(numberOfGausspoints), ')!']);
            % sum of all first derivatives must equal 0
            dNr = reshape(dN_xi_k_I, [], size(dN_xi_k_I, 3));
            assert(all(sum(dNr, 2)-ones(dimension*numberOfGausspoints, 1) < tolerance), ['Error in first derivatives of Lagrange shape function (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ', numberOfGausspoints=', num2str(numberOfGausspoints), ')!']);
            % sum of all second derivatives must equal 0
            d2Nr = reshape(d2N_xi_xi_k_I, [], size(d2N_xi_xi_k_I, 4));
            assert(all(sum(d2Nr, 2)-ones(dimension^2*numberOfGausspoints, 1) < tolerance), ['Error in second derivatives of Lagrange shape function (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ', numberOfGausspoints=', num2str(numberOfGausspoints), ')!']);
        end
    end
end
disp('Symbolic computation of Lagrange shape functions completed successfully.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validate Lagrange shape functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for dimension = 1:3
    numberOfNodesArray = numberOfNodesCell{dimension};
    numberOfGausspointsArray = numberOfGausspointsCell{dimension};
    for ii = 1:length(numberOfNodesArray)
        numberOfNodes = numberOfNodesArray(ii);
        elementGeometryType = determineElementGeometryType(dimension, numberOfNodes);
        for jj = 1:length(numberOfGausspointsArray)
            numberOfGausspoints = numberOfGausspointsArray(jj);
            gaussPoints = gaussPointsAndWeights(dimension, numberOfGausspoints, elementGeometryType);
            [N_k_I, dN_xi_k_I, d2N_xi_xi_k_I] = computeLagrangeShapeFunction(dimension, numberOfNodes, numberOfGausspoints, gaussPoints);
            NSymbolic_k_I = lagrangeShapeFunctionCell{dimension, ii, jj}.N_k_I;
            dNSymbolic_xi_k_I = lagrangeShapeFunctionCell{dimension, ii, jj}.dN_xi_k_I;
            d2NSymbolic_xi_xi_k_I = lagrangeShapeFunctionCell{dimension, ii, jj}.d2N_xi_xi_k_I;

            assert(max(max(abs(N_k_I-NSymbolic_k_I))) < tolerance, ['Error in Lagrange shape function (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ', numberOfGausspoints=', num2str(numberOfGausspoints), ')!']);
            assert(max(max(max(abs(dN_xi_k_I-dNSymbolic_xi_k_I)))) < tolerance, ['Error in first derivatives of Lagrange shape function (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ', numberOfGausspoints=', num2str(numberOfGausspoints), ')!']);
            assert(max(max(max(max(abs(d2N_xi_xi_k_I-d2NSymbolic_xi_xi_k_I))))) < tolerance, ['Error in second derivatives of Lagrange shape function (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ', numberOfGausspoints=', num2str(numberOfGausspoints), ')!']);
        end
    end
end
disp('Validation of Lagrange shape functions completed successfully.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validate metric shape functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for dimension = 1:3
    numberOfNodesArray = numberOfNodesMetricShapeFunctionsCell{dimension};
    numberOfGausspointsArray = numberOfGausspointsCell{dimension};
    for ii = 1:length(numberOfNodesArray)
        numberOfNodes = numberOfNodesArray(ii);
        continuumObject.elementGeometryType = determineElementGeometryType(dimension, numberOfNodes);
        nodalPoints = elementNodesInLocalCoordinates(dimension, continuumObject.elementGeometryType, numberOfNodes);
        indexNumberOfNodes = find(numberOfNodesCell{dimension} == numberOfNodes);
        for jj = 1:length(numberOfGausspointsArray)
            numberOfGausspoints = numberOfGausspointsArray(jj);
            gaussPoints = gaussPointsAndWeights(dimension, numberOfGausspoints, continuumObject.elementGeometryType);
            [M_k_I, dMr] = computeMetricShapeFunctions(continuumObject, dimension, nodalPoints, gaussPoints);
            N_k_I = lagrangeShapeFunctionCell{dimension, indexNumberOfNodes, jj}.N_k_I;
            dN_xi_k_I = lagrangeShapeFunctionCell{dimension, indexNumberOfNodes, jj}.dN_xi_k_I;
            dNr = reshape(dN_xi_k_I, [], size(dN_xi_k_I, 3));

            assert(max(max(abs(M_k_I-N_k_I))) < tolerance, ['Error in metric shape function (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ', numberOfGausspoints=', num2str(numberOfGausspoints), ')!']);
            assert(max(max(abs(dMr-dNr))) < tolerance, ['Error in first derivatives of metric shape function (dimension=', num2str(dimension), ', numberOfNodes=', num2str(numberOfNodes), ', numberOfGausspoints=', num2str(numberOfGausspoints), ')!']);
        end
    end
end
disp('Validation of metric shape functions completed successfully.');
