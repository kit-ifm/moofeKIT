function [M_k_I, dMr] = computeMetricShapeFunctions(continuumObject, dimension, nodalPoints, evaluationPoints)
%COMPUTEMETRICSHAPEFUNCTIONS Computes metric shape functions.
%   This function computes the metric shape functions commonly used in the
%   context of Petrov-Galerkin Finite Elements
%
%   CALL
%   [M, dMr] = computeMetricShapeFunctions(continuumObject, dimension, nodalPoints, evaluationPoints)
%   continuumObject: object of solidSuperClass
%   dimension: dimension
%   nodalPoints: nodal points in physical (x, y) or skew coordinates (xBar, yBar)
%   evaluationPoints: evaluation points (e.g. gauss points) in physical or skew coordinates
%   M_k_I: Shape functions
%   dMr:  Derivatives of the shape functions with respect to the given physical or skew coordinates
%
%   REFERENCES
%   https://doi.org/10.1007/s10999-014-9289-3
%   https://doi.org/10.1002/nme.6817
%
%   CREATOR(S)
%   Felix Zaehringer

% check & normalize input
assert(dimension <= 3, 'Not implemted yet for given dimension!');
assert(size(nodalPoints, 1) == dimension || size(nodalPoints, 2) == dimension, 'Size of nodalPoints matrix does not fit the given dimension!');
assert(size(evaluationPoints, 1) == dimension || size(evaluationPoints, 2) == dimension, 'Size of evaluationPoints matrix does not fit the given dimension!');
if size(nodalPoints, 2) == dimension
    nodalPoints = nodalPoints.';
end
if size(evaluationPoints, 2) == dimension
    evaluationPoints = evaluationPoints.';
end

% compute metric shape functions
numberOfNodes = size(nodalPoints, 2);
P = computePMatrix(continuumObject, nodalPoints, numberOfNodes, dimension);
p = computePMatrix(continuumObject, evaluationPoints, numberOfNodes, dimension);
M_k_I = p / P;

% compute derivatives of metric shape functions
q = computeQMatrix(continuumObject, evaluationPoints, numberOfNodes, dimension);
dMr = q / P;

end

function P = computePMatrix(continuumObject, points, numberOfNodes, dimension)
% P Matrix in analogy to https://doi.org/10.1007/s10999-014-9289-3
numberOfPoints = size(points, 2);
xPoints = points(1, :).';
if dimension >= 2
    yPoints = points(2, :).';
end
if dimension >= 3
    zPoints = points(3, :).';
end
onesColumn = ones(numberOfPoints, 1);
if dimension == 1
    if strcmp(continuumObject.elementGeometryType, 'oneDimensional')
        if numberOfNodes == 2
            P = [onesColumn, xPoints];
        elseif numberOfNodes == 3
            P = [onesColumn, xPoints, xPoints.^2];
        else
            error(['Not implemented yet for given numberOfNodes (dimension=1, numberOfNodes=', num2str(numberOfNodes),')!']);
        end
    else
        error('Not implemented yet for given elementGeometryType!');
    end
elseif dimension == 2
    if strcmp(continuumObject.elementGeometryType, 'quadrilateral')
        if numberOfNodes == 4
            P = [onesColumn, xPoints, yPoints, xPoints .* yPoints];
        elseif numberOfNodes == 8
            P = [onesColumn, xPoints, yPoints, xPoints.^2, xPoints .* yPoints, yPoints.^2, xPoints.^2 .* yPoints, xPoints .* yPoints.^2];
        elseif numberOfNodes == 9
            P = [onesColumn, xPoints, yPoints, xPoints.^2, xPoints .* yPoints, yPoints.^2, xPoints.^2 .* yPoints, xPoints .* yPoints.^2, xPoints.^2 .* yPoints.^2];
        else
            error(['Not implemented yet for given numberOfNodes (dimension=2, numberOfNodes=', num2str(numberOfNodes),')!']);
        end
    else
        error('Not implemented yet for given elementGeometryType!');
    end
elseif dimension == 3
    if strcmp(continuumObject.elementGeometryType, 'hexahedral')
        if numberOfNodes == 8
            P = [onesColumn, xPoints, yPoints, zPoints, yPoints .* zPoints, zPoints .* xPoints, xPoints .* yPoints, xPoints .* yPoints .* zPoints];
        elseif numberOfNodes == 20
            P = [onesColumn, xPoints, yPoints, zPoints, xPoints.^2, yPoints.^2, zPoints.^2, yPoints .* zPoints, zPoints .* xPoints, xPoints .* yPoints, xPoints.^2 .* yPoints, xPoints.* yPoints.^2, yPoints.^2 .* zPoints, yPoints .* zPoints.^2, zPoints.^2 .* xPoints, xPoints.^2 .* zPoints, xPoints .* yPoints .* zPoints, xPoints.^2 .* yPoints .* zPoints, xPoints .* yPoints.^2 .* zPoints, xPoints .* yPoints .* zPoints.^2];
        else
            error(['Not implemented yet for given numberOfNodes (dimension=3, numberOfNodes=', num2str(numberOfNodes),')!']);
        end
    else
        error('Not implemented yet for given elementGeometryType!');
    end
else
    error('Not implemented yet for given dimension!');
end
end

function Q = computeQMatrix(continuumObject, evaluationPoints, numberOfNodes, dimension)
% Q Matrix => derivative of the P Matrix
numberOfEvaluationPoints = size(evaluationPoints, 2);
Q = zeros(dimension*numberOfEvaluationPoints, numberOfNodes);
xPoints = evaluationPoints(1, :).';
if dimension >= 2
    yPoints = evaluationPoints(2, :).';
end
if dimension >= 3
    zPoints = evaluationPoints(3, :).';
end
if dimension == 1
    if strcmp(continuumObject.elementGeometryType, 'oneDimensional')
        if numberOfNodes == 2
            Q(1:dimension:end, 2) = 1;
        elseif numberOfNodes == 3
            Q(1:dimension:end, 2) = 1;
            Q(1:dimension:end, 3) = 2 * xPoints;
        else
            error('Not implemented yet for given numberOfNodes!');
        end
    else
        error('Not implemented yet for given elementGeometryType!');
    end
elseif dimension == 2
    if strcmp(continuumObject.elementGeometryType, 'quadrilateral')
        if numberOfNodes == 4
            Q(1:dimension:end, 2) = 1;
            Q(2:dimension:end, 3) = 1;
            Q(1:dimension:end, 4) = yPoints;
            Q(2:dimension:end, 4) = xPoints;
        elseif numberOfNodes == 8 || numberOfNodes == 9
            Q(1:dimension:end, 2) = 1;
            Q(2:dimension:end, 3) = 1;
            Q(1:dimension:end, 4) = 2 * xPoints;
            Q(1:dimension:end, 5) = yPoints;
            Q(2:dimension:end, 5) = xPoints;
            Q(2:dimension:end, 6) = 2 * yPoints;
            Q(1:dimension:end, 7) = 2 * xPoints .* yPoints;
            Q(2:dimension:end, 7) = xPoints.^2;
            Q(1:dimension:end, 8) = yPoints.^2;
            Q(2:dimension:end, 8) = 2 * xPoints .* yPoints;
            if numberOfNodes == 9
                Q(1:dimension:end, 9) = 2 * xPoints .* yPoints.^2;
                Q(2:dimension:end, 9) = 2 * xPoints.^2 .* yPoints;
            end
        else
            error('Not implemented yet for given numberOfNodes!');
        end
    else
        error('Not implemented yet for given elementGeometryType!');
    end
elseif dimension == 3
    if strcmp(continuumObject.elementGeometryType, 'hexahedral')
        if numberOfNodes == 8
            Q(1:dimension:end, 2) = 1;
            Q(2:dimension:end, 3) = 1;
            Q(3:dimension:end, 4) = 1;
            Q(2:dimension:end, 5) = zPoints;
            Q(3:dimension:end, 5) = yPoints;
            Q(1:dimension:end, 6) = zPoints;
            Q(3:dimension:end, 6) = xPoints;
            Q(1:dimension:end, 7) = yPoints;
            Q(2:dimension:end, 7) = xPoints;
            Q(1:dimension:end, 8) = yPoints .* zPoints;
            Q(2:dimension:end, 8) = xPoints .* zPoints;
            Q(3:dimension:end, 8) = xPoints .* yPoints;
        elseif numberOfNodes == 20
            Q(1:dimension:end, 2) = 1;
            Q(2:dimension:end, 3) = 1;
            Q(3:dimension:end, 4) = 1;
            Q(1:dimension:end, 5) = 2 * xPoints;
            Q(2:dimension:end, 6) = 2 * yPoints;
            Q(3:dimension:end, 7) = 2 * zPoints;
            Q(2:dimension:end, 8) = zPoints;
            Q(3:dimension:end, 8) = yPoints;
            Q(1:dimension:end, 9) = zPoints;
            Q(3:dimension:end, 9) = xPoints;
            Q(1:dimension:end, 10) = yPoints;
            Q(2:dimension:end, 10) = xPoints;
            Q(1:dimension:end, 11) = 2 * xPoints .* yPoints;
            Q(2:dimension:end, 11) = xPoints.^2;
            Q(1:dimension:end, 12) = yPoints.^2;
            Q(2:dimension:end, 12) = 2 * xPoints .* yPoints;
            Q(2:dimension:end, 13) = 2 * yPoints .* zPoints;
            Q(3:dimension:end, 13) = yPoints.^2;
            Q(2:dimension:end, 14) = zPoints.^2;
            Q(3:dimension:end, 14) = 2 * yPoints .* zPoints;
            Q(1:dimension:end, 15) = zPoints.^2;
            Q(3:dimension:end, 15) = 2 * xPoints .* zPoints;
            Q(1:dimension:end, 16) = 2 * xPoints .* zPoints;
            Q(3:dimension:end, 16) = xPoints.^2;
            Q(1:dimension:end, 17) = yPoints .* zPoints;
            Q(2:dimension:end, 17) = xPoints .* zPoints;
            Q(3:dimension:end, 17) = xPoints .* yPoints;
            Q(1:dimension:end, 18) = 2 * xPoints .* yPoints .* zPoints;
            Q(2:dimension:end, 18) = xPoints.^2 .* zPoints;
            Q(3:dimension:end, 18) = xPoints.^2 .* yPoints;
            Q(1:dimension:end, 19) = yPoints.^2 .* zPoints;
            Q(2:dimension:end, 19) = 2 * xPoints .* yPoints .* zPoints;
            Q(3:dimension:end, 19) = yPoints.^2 .* xPoints;
            Q(1:dimension:end, 20) = yPoints .* zPoints.^2;
            Q(2:dimension:end, 20) = xPoints .* zPoints.^2;
            Q(3:dimension:end, 20) = 2 * xPoints .* yPoints .* zPoints;
        else
            error('Not implemented yet for given numberOfNodes!');
        end
    else
        error('Not implemented yet for given elementGeometryType!');
    end
else
    error('Not implemented yet for given dimension!');
end
end
