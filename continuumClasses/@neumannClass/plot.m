function plot(obj, dofObject, setupObject)
%plots loads of neumann
plotObject = setupObject.plotObject;
time = plotObject.time;

%formatting
lineWeight = plotObject.lineWeight;
deltaXY = plotObject.deltaXY;
lengthPoint = 0.2;
phiPoint = pi / 6;

%% getting nodal data
updateNodalForces(obj, dofObject);
nodalForces = obj.nodalForces;

maxF = max(max(abs(nodalForces)));
switch time
    case 'R'
        nodes = obj.masterObject.qR;
    case 'N'
        nodes = obj.masterObject.qN;
    case 'N1'
        nodes = obj.masterObject.qN1;
    otherwise
        error('unknown timepoint')
end
if isa(obj.masterObject, 'plateClass')
    nodes = [obj.masterObject.meshObject.nodes, zeros(size(obj.masterObject.meshObject.nodes, 1), 1)];
    switch upper(time)
        case 'R'
            nodes(:, 3) = nodes(:, 3) + obj.masterObject.qR(:, 1);
        case 'N'
            nodes(:, 3) = nodes(:, 3) + obj.masterObject.qN(:, 1);
        case 'N1'
            nodes(:, 3) = nodes(:, 3) + obj.masterObject.qN1(:, 1);
        otherwise
            error('unknown timepoint')
    end
end

%adjust dimension of arrays if 1D: everything in y-direction is zero
if obj.masterObject.dimension == 1
    nodalForces = [nodalForces, zeros(size(nodalForces, 1), 1)];
    nodes = [nodes, zeros(size(nodes, 1), 1)];
end

numLoad = size(nodalForces, 1);

%% plotting of mesh
if ~isa(obj.masterObject, 'plateClass')
    % 2D solids
    x = zeros(numLoad, 5);
    y = zeros(numLoad, 5);

    for e = 1:numLoad
        Fx = nodalForces(e, 1);
        Fy = nodalForces(e, 2);
        F = sqrt(Fx^2+Fy^2);

        if F > eps
            points = zeros(2, 5);
            L = F / maxF * 0.08 * deltaXY;

            points(1, 1) = L;
            points(1, 2) = 0;
            points(1, 3) = L * lengthPoint;
            points(1, 4) = 0;
            points(1, 5) = L * lengthPoint;

            points(2, 1) = 0;
            points(2, 2) = 0;
            points(2, 3) = L * lengthPoint * tan(phiPoint);
            points(2, 4) = 0;
            points(2, 5) = -L * lengthPoint * tan(phiPoint);

            %rotation and translation
            sinP = Fy / F;
            cosP = -Fx / F;
            A = [cosP, sinP; -sinP, cosP];
            points = A * points;

            x(e, :) = points(1, :) + nodes(obj.masterNodes(e), 1);
            y(e, :) = points(2, :) + nodes(obj.masterNodes(e), 2);
        end
    end
    plot(x', y', 'b', 'linewidth', lineWeight)
else
    % plates
    x = zeros(numLoad, 1);
    y = zeros(numLoad, 1);
    z = zeros(numLoad, 1);
    w = zeros(numLoad, 1);
    for e = 1:numLoad
        Fz = nodalForces(e, 1);

        if abs(Fz) > eps
            L = Fz / maxF * 0.08 * deltaXY;
            x(e) = nodes(obj.masterNodes(e), 1);
            y(e) = nodes(obj.masterNodes(e), 2);
            z(e) = nodes(obj.masterNodes(e), 3) - L;
            w(e) = L;
        end
    end
    quiver3(x, y, z, zeros(size(x)), zeros(size(x)), w, 'off', 'b', 'linewidth', lineWeight);
end
end