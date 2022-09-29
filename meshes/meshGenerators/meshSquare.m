function [NODES,EDOF,bounEDOFs] = meshSquare(numberOfElements, order, x0, xEnd)
    [NODES,EDOF,bounEDOFs] = meshRectangle(xEnd-x0, xEnd-x0, numberOfElements, numberOfElements, order);
    NODES = NODES+(xEnd-x0)/2+x0;
end