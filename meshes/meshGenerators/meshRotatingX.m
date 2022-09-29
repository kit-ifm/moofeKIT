function [nodes, edof] = meshRotatingX(armLength, armWidth, bodyThickness, numberOfElementsArmLength, numberOfElementsArmWidth, numberOfElementsThickness, order, serendipity)
%MESHROTATINGX Returns the mesh of an X-shaped body
%
% CALL
% meshRotatingX(bodyThickness, armLength, armWidth, numberOfElementsBodyThickness, numberOfElementsArmLength, numberOfElementsArmWidth)
% bodyThickness: thickness of the X-shaped body
% armLength: length of each arm
% armWidth: width of each arm
% numberOfElementsThickness: number of finite elements in the direction
% of the body thickness.
% numberOfElementsArmLength: number of finite elements in length direction
% of each arm
% numberOfElementsArmWidth: number of finite elements in width direction
% of each arm
% order: interpolation order of the elements (may be 1 or 2)
% serendipity: if serendipity elements are to be used for order==2
%
% REFERENCE
% In preparation
%
% CREATOR(S)
% Felix Zaehringer

% set standard arguments
if nargin <= 6
    order = 1;
end
if nargin <= 7
    serendipity = false;
end
if order == 1 && serendipity == true
    warning('setting serendipity=true has no effect if order is 1')
end

% generate mesh
tolerance = 1e-9;
[middlePartNodes, middlePartEdof] = meshGeneratorCube(armWidth, armWidth, bodyThickness, numberOfElementsArmWidth, numberOfElementsArmWidth, numberOfElementsThickness, order, serendipity);
[leftArmNodes, leftArmEdof] = meshGeneratorCube(armLength, armWidth, bodyThickness, numberOfElementsArmLength, numberOfElementsArmWidth, numberOfElementsThickness, order, serendipity);
[rightArmNodes, rightArmEdof] = meshGeneratorCube(armLength, armWidth, bodyThickness, numberOfElementsArmLength, numberOfElementsArmWidth, numberOfElementsThickness, order, serendipity);
[topArmNodes, topArmEdof] = meshGeneratorCube(armWidth, armLength, bodyThickness, numberOfElementsArmWidth, numberOfElementsArmLength, numberOfElementsThickness, order, serendipity);
[bottomArmNodes, bottomArmEdof] = meshGeneratorCube(armWidth, armLength, bodyThickness, numberOfElementsArmWidth, numberOfElementsArmLength, numberOfElementsThickness, order, serendipity);
leftArmNodes(:, 1) = leftArmNodes(:, 1) - (armWidth + armLength) / 2;
rightArmNodes(:, 1) = rightArmNodes(:, 1) + (armWidth + armLength) / 2;
topArmNodes(:, 2) = topArmNodes(:, 2) + (armWidth + armLength) / 2;
bottomArmNodes(:, 2) = bottomArmNodes(:, 2) - (armWidth + armLength) / 2;

leftArmEdof = leftArmEdof + max(max(middlePartEdof));
[nodes, edof] = meshTie(middlePartNodes, middlePartEdof, leftArmNodes, leftArmEdof, tolerance);
rightArmEdof = rightArmEdof + max(max(edof));
[nodes, edof] = meshTie(nodes, edof, rightArmNodes, rightArmEdof, tolerance);
topArmEdof = topArmEdof + max(max(edof));
[nodes, edof] = meshTie(nodes, edof, topArmNodes, topArmEdof, tolerance);
bottomArmEdof = bottomArmEdof + max(max(edof));
[nodes, edof] = meshTie(nodes, edof, bottomArmNodes, bottomArmEdof, tolerance);
end
