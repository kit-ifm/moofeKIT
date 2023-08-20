function [nodes, edof, bounEdof] = meshPatchTestDistorted2D(lengthX, lengthY, order, serendipity)
%MESHPATCHTESTDISTORTED Initially distorted mesh used in the patch test
%   This function returns the nodes and the edof for the initially
%   distorted mesh that is used in the patch test.
%
%   REFERENCE
%   https://doi.org/10.1016/0168-874X(85)90003-4
%
%   AUTHOR(S)
%   Felix Zaehringer

if order == 1 && ~serendipity
    nodes = [0, 0.24, 0.24, 0, 0.04, 0.18, 0.16, 0.08; ...
        0, 0, 0.12, 0.12, 0.02, 0.03, 0.08, 0.08]';
    nodes(:, 1) = nodes(:, 1) * lengthX / 0.24;
    nodes(:, 2) = nodes(:, 2) * lengthY / 0.12;
    edof = [1, 2, 6, 5; ...
        2, 3, 7, 6; ...
        8, 7, 3, 4; ...
        1, 5, 8, 4; ...
        5, 6, 7, 8];
    bounEdof = struct('SX1', [4, 1], 'SX2', [2, 3], 'SY1', [1, 2], 'SY2', [3, 4]);
elseif order == 2
    if serendipity
        nodes = [0, 0.24, 0.24, 0, 0.04, 0.18, 0.16, 0.08, 0.12, 0.24, 0.12, 0, 0.02, 0.11, 0.21, 0.17, 0.2, 0.12, 0.04, 0.06; ...
            0, 0, 0.12, 0.12, 0.02, 0.03, 0.08, 0.08, 0, 0.06, 0.12, 0.06, 0.01, 0.025, 0.015, 0.055, 0.1, 0.08, 0.1, 0.05]';
        edof = [1, 2, 6, 5, 9, 15, 14, 13; ...
            2, 3, 7, 6, 10, 17, 16, 15; ...
            8, 7, 3, 4, 18, 17, 11, 19; ...
            1, 5, 8, 4, 13, 20, 19, 12; ...
            5, 6, 7, 8, 14, 16, 18, 20];
    else
        nodes = [0, 0.24, 0.24, 0, 0.04, 0.18, 0.16, 0.08, 0.12, 0.24, 0.12, 0, 0.02, 0.11, 0.21, 0.17, 0.2, 0.12, 0.04, 0.06, 0.115, 0.205, 0.12, 0.03, 0.115; ...
            0, 0, 0.12, 0.12, 0.02, 0.03, 0.08, 0.08, 0, 0.06, 0.12, 0.06, 0.01, 0.025, 0.015, 0.055, 0.1, 0.08, 0.1, 0.05, 0.0125, 0.0575, 0.1, 0.055, 0.0525]';
        edof = [1, 2, 6, 5, 9, 15, 14, 13, 21; ...
            2, 3, 7, 6, 10, 17, 16, 15, 22; ...
            8, 7, 3, 4, 18, 17, 11, 19, 23; ...
            1, 5, 8, 4, 13, 20, 19, 12, 24; ...
            5, 6, 7, 8, 14, 16, 18, 20, 25];
    end
    nodes(:, 1) = nodes(:, 1) * lengthX / 0.24;
    nodes(:, 2) = nodes(:, 2) * lengthY / 0.12;
    bounEdof = struct('SX1', [4, 12, 1], 'SX2', [2, 10, 3], 'SY1', [1, 9, 2], 'SY2', [3, 11, 4]);
else
    error('Not implemented!')
end
end
