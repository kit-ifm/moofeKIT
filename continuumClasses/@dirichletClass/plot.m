function plot(obj, dofObject, setupObject)
%plots boundaries from dirichlet

plotObject = setupObject.plotObject;
time = plotObject.time;

nodeListDirichletBoundary = obj.nodeList;
numberOfNodesDirichletBoundary = size(nodeListDirichletBoundary, 1);

%formatting
a = plotObject.boundarySize * plotObject.deltaXY;
h = a * sqrt(3) / 2;
lineWeight = plotObject.lineWeight;

%% getting nodal data
switch upper(time)
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
    nodes = [nodes, zeros(size(nodes, 1), 1)];
end

%% plotting
x = zeros(numberOfNodesDirichletBoundary, 4);
y = zeros(numberOfNodesDirichletBoundary, 4);
z = zeros(numberOfNodesDirichletBoundary, 4);

if isa(obj.masterObject, 'plateClass')
    % plates
    for i = 1:numel(obj.nodalDof)
        switch obj.nodalDof(i)
            case 1
                x1 = zeros(numberOfNodesDirichletBoundary, 4);
                y1 = zeros(numberOfNodesDirichletBoundary, 4);
                x2 = zeros(numberOfNodesDirichletBoundary, 4);
                y2 = zeros(numberOfNodesDirichletBoundary, 4);
                z = zeros(numberOfNodesDirichletBoundary, 4);

                xyTri = [0, a / 2, -a / 2, 0];
                zTri = [0, -h, -h, 0];

                x1(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), xyTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 1));
                y1(:, :) = kron(ones(1, 4), nodes(nodeListDirichletBoundary, 2));
                x2(:, :) = kron(ones(1, 4), nodes(nodeListDirichletBoundary, 1));
                y2(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), xyTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 2));
                z(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), zTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 3));

                plot3(x1', y1', z', 'r', 'linewidth', 2*lineWeight)
                plot3(x2', y2', z', 'r', 'linewidth', 2*lineWeight)
            case 2
                x = zeros(numberOfNodesDirichletBoundary, 5);
                y = zeros(numberOfNodesDirichletBoundary, 5);
                z = zeros(numberOfNodesDirichletBoundary, 5);

                xyTri = [-a / 2, a / 2, a / 2, -a / 2, -a / 2];
                zTri = [-a / 2, -a / 2, a / 2, a / 2, -a / 2];

                x(:, :) = kron(ones(1, 5), nodes(nodeListDirichletBoundary, 1));
                y(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), xyTri) + kron(ones(1, 5), nodes(nodeListDirichletBoundary, 2));
                z(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), zTri) + kron(ones(1, 5), nodes(nodeListDirichletBoundary, 3));

                plot3(x', y', z', 'r', 'linewidth', 2*lineWeight)
            case 3
                x = zeros(numberOfNodesDirichletBoundary, 5);
                y = zeros(numberOfNodesDirichletBoundary, 5);
                z = zeros(numberOfNodesDirichletBoundary, 5);

                xyTri = [-a / 2, a / 2, a / 2, -a / 2, -a / 2];
                zTri = [-a / 2, -a / 2, a / 2, a / 2, -a / 2];

                x(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), xyTri) + kron(ones(1, 5), nodes(nodeListDirichletBoundary, 1));
                y(:, :) = kron(ones(1, 5), nodes(nodeListDirichletBoundary, 2));
                z(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), zTri) + kron(ones(1, 5), nodes(nodeListDirichletBoundary, 3));

                plot3(x', y', z', 'r', 'linewidth', 2*lineWeight)
            otherwise
                error('direction not implemented')
        end
    end
elseif obj.masterObject.dimension == 2
    % general 2D solids
    for i = 1:numel(obj.nodalDof)
        switch obj.nodalDof(i)
            case 1
                xTri = [0, -h, -h, 0];
                yTri = [0, a / 2, -a / 2, 0];
            case 2
                xTri = [0, -a / 2, a / 2, 0];
                yTri = [0, -h, -h, 0];
            otherwise
                if any(strcmp(class(obj.masterObject), {'solidThermoClass', 'solidVelocityClass'}))
                    break;
                end
                error('direction not implemented')
        end

        x(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), xTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 1));
        y(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), yTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 2));
        plot(x', y', 'r', 'linewidth', 2*lineWeight)
    end
elseif obj.masterObject.dimension == 3 && any(strcmp(class(obj.masterObject), {'solidClass', 'solidVelocityClass', 'solidShellClass'}))
    % 3D solids
    for i = 1:numel(obj.nodalDof)
        switch obj.nodalDof(i)
            case 1
                xTri = [0, -h, -h, 0];
                yTri = [0, a / 2, -a / 2, 0];
                zTri = [0, a / 2, -a / 2, 0];

                x(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), xTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 1));
                y(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), yTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 2));
                z(:, :) = kron(ones(1, 4), nodes(nodeListDirichletBoundary, 3));
                plot3(x', y', z', 'r', 'linewidth', 2*lineWeight)

                x(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), xTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 1));
                y(:, :) = kron(ones(1, 4), nodes(nodeListDirichletBoundary, 2));
                z(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), zTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 3));
                plot3(x', y', z', 'r', 'linewidth', 2*lineWeight)

            case 2
                xTri = [0, -a / 2, a / 2, 0];
                yTri = [0, -h, -h, 0];
                zTri = [0, a / 2, -a / 2, 0];
                    
                x(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), xTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 1));
                y(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), yTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 2));
                z(:, :) = kron(ones(1, 4), nodes(nodeListDirichletBoundary, 3));
                plot3(x', y', z', 'r', 'linewidth', 2*lineWeight)

                x(:, :) = kron(ones(1, 4), nodes(nodeListDirichletBoundary, 1));
                y(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), yTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 2));
                z(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), zTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 3));
                plot3(x', y', z', 'r', 'linewidth', 2*lineWeight)

            case 3
                xTri = [0, a / 2, -a / 2, 0];
                yTri = [0, a / 2, -a / 2, 0];
                zTri = [0, -h, -h, 0];

                x(:, :) = kron(ones(1, 4), nodes(nodeListDirichletBoundary, 1));
                y(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), yTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 2));
                z(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), zTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 3));
                plot3(x', y', z', 'r', 'linewidth', 2*lineWeight)
    
                x(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), xTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 1));
                y(:, :) = kron(ones(1, 4), nodes(nodeListDirichletBoundary, 2));
                z(:, :) = kron(ones(numberOfNodesDirichletBoundary, 1), zTri) + kron(ones(1, 4), nodes(nodeListDirichletBoundary, 3));
                plot3(x', y', z', 'r', 'linewidth', 2*lineWeight)

                otherwise
                % error('direction not implemented')
        end    
    end
end