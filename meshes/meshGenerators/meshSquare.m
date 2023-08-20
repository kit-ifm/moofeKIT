function [NODES, EDOF, bounEDOFs] = meshSquare(numberOfElements, order, x0, xEnd, serendipity)
%set standard arguments
if nargin <= 4
    serendipity = false;
end

[NODES, EDOF, bounEDOFs] = meshRectangle(xEnd-x0, xEnd-x0, numberOfElements, numberOfElements, order, serendipity);
NODES = NODES + (xEnd - x0) / 2 + x0;
end