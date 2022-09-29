function dataFE = massMatrix(obj)
dataFE = [];
for index = 1:numel(obj)
    objIndex = obj(index);
    globalFullEdof = objIndex.globalFullEdof;
    edof = objIndex.edof; 
    numberOfGausspoints = objIndex.numberOfGausspoints;
    gaussWeight = objIndex.shapeFunctions.gaussWeight;
    NAll = objIndex.shapeFunctions.N';
    dNrAll = objIndex.shapeFunctions.dNr';
    rho = objIndex.materialData.rho;
    QR = objIndex.QR;
    numberOfElements = size(globalFullEdof,1);
    dimension = objIndex.dimension;
    dataFEindex(numberOfElements) = struct('MeVector',[],'indexMi',[],'indexMj',[]);
    parfor e = 1:numberOfElements
        [globalEdofH1,globalEdofH2] = expandEdof(globalFullEdof(e,:));
        dataFEindex(e).indexMi = globalEdofH1';
        dataFEindex(e).indexMj = globalEdofH2';
        numberOfDofs = numel(globalFullEdof(e,:));        
        dNr = dNrAll;
        N = NAll;
        sizeN = size(N,2)*(dimension);
        Me = zeros(numberOfDofs);
        J = QR(edof(e,:),1:dimension)'*dNr;
        for gp = 1:numberOfGausspoints
            index2 = dimension*gp-(dimension-1):dimension*gp;
            detJ = det(J(:,index2)');
            if detJ <= 0
                error('Jacobideterminant less or equal zero.')
            end
            A1 = (N(gp,:)'*N(gp,:))*rho;
            MAT = zeros(numberOfDofs);
            for k = 1:dimension
                MAT(k:dimension:sizeN,k:dimension:sizeN) = A1;
            end
            Me = Me + MAT*detJ*gaussWeight(gp);
        end
        dataFEindex(e).MeVector = Me(:);
    end
    dataFE = [dataFE dataFEindex];
end