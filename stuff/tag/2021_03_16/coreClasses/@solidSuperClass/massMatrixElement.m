function dataFE = massMatrixElement(obj)
globalFullEdof = obj.globalFullEdof;
edof = obj.edof;
numberOfGausspoints = obj.numberOfGausspoints;
gaussWeight = obj.shapeFunctions.gaussWeight;
NAll = obj.shapeFunctions.N';
dNrAll = obj.shapeFunctions.dNr';
rho = obj.materialData.rho;
qR = obj.qR;
numberOfElements = size(globalFullEdof,1);
dimension = obj.dimension;
additionalFields = obj.additionalFields;
dataFE(numberOfElements) = struct('MeVector',[],'indexMi',[],'indexMj',[]);
parfor e = 1:numberOfElements
    [globalEdofH1,globalEdofH2] = expandEdof(globalFullEdof(e,:));
    dataFE(e).indexMi = globalEdofH1';
    dataFE(e).indexMj = globalEdofH2';
    numberOfDofs = numel(globalFullEdof(e,:));
    dNr = dNrAll;
    N = NAll;    
    sizeN = size(N,2)*(dimension + additionalFields);
    Me = zeros(numberOfDofs);
    J = qR(edof(e,:),1:dimension)'*dNr;
    for gp = 1:numberOfGausspoints
        index2 = dimension*gp-(dimension-1):dimension*gp;
        detJ = det(J(:,index2)');
        if detJ <= 0
            error('Jacobideterminant less or equal zero.')
        end
        A1 = (N(gp,:)'*N(gp,:))*rho;
        MAT = zeros(numberOfDofs);
        for k = 1:dimension
            MAT(k:dimension+additionalFields:sizeN,k:dimension+additionalFields:sizeN) = A1;
        end
        Me = Me + MAT*detJ*gaussWeight(gp);
    end
    dataFE(e).MeVector = Me(:);
end
