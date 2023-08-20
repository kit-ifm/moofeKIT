classdef storageFEClass < matlab.mixin.Copyable
    properties
        dataFE % only for neumann boundary
    end
    methods
        function [dataFE, initialDataFE] = initializeDataFE(obj, requiredData)
            if strcmp(requiredData, 'residualAndTangent')
                dataFE = struct('indexReI', [], 'Re', [], 'indexKeI', [], 'indexKeJ', [], 'Ke', []);
            elseif strcmp(requiredData, 'postData')
                dataFE = struct('indexReI', [], 'Se', [], 'indexKeI', [], 'indexKeJ', [], 'Me', []);
            elseif strcmp(requiredData, 'massMatrix')
                dataFE = struct('Me', [], 'indexMeI', [], 'indexMeJ', []);
            end
            initialDataFE = dataFE;
        end
        function dataFEMixed = initializeDataFEMixed(obj)
            dataFEMixed = struct('ReDTilde', [], 'KeDDTilde', []);
        end
        function [array, stressTensor, arrayMixed] = initializeArrayStress(obj, numberOfDofs, numberOfNodes)
            array = struct('Re', zeros(numberOfDofs, 1), 'Ke', zeros(numberOfDofs), ... % solver: computing residual and tangent
                'Se', zeros(numberOfNodes, 1), 'Me', zeros(numberOfNodes)); % postprocessing: computing stress etc.
            stressTensor = struct('FirstPK', 'Cauchy');
        end
        function dataFE = assignElementDataFE(obj, dataFE, array, globalFullEdof, e, dimension, additionalFields, N, requiredData)
            if strcmp(requiredData, 'residualAndTangent')
                dataFE.indexReI = double(globalFullEdof(e, :))';
                [globalFullEdof1, globalFullEdof2] = expandEdof(globalFullEdof(e, :));
                dataFE.indexKeI = globalFullEdof1';
                dataFE.indexKeJ = globalFullEdof2';
                dataFE.Re = array.Re;
                dataFE.Ke = array.Ke(:);
            elseif strcmp(requiredData, 'postData')
                globalFullEdof = globalFullEdof(:, 1:dimension+additionalFields:(dimension + additionalFields)*size(N, 2));
                dataFE.indexReI = double(globalFullEdof(e, :))';
                [globalFullEdof1, globalFullEdof2] = expandEdof(globalFullEdof(e, :));
                dataFE.indexKeI = globalFullEdof1';
                dataFE.indexKeJ = globalFullEdof2';
                dataFE.Se = array.Se;
                dataFE.Me = array.Me(:);
            elseif strcmp(requiredData, 'massMatrix')
                [globalFullEdof1, globalFullEdof2] = expandEdof(globalFullEdof(e, :));
                dataFE.indexMeI = globalFullEdof1';
                dataFE.indexMeJ = globalFullEdof2';
                dataFE.Me = array.Me(:);
            end
        end
    end
end
