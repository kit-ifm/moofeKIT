classdef storageFEClass < matlab.mixin.Copyable
    properties
        dataFE % only for neumann boundary
    end
    methods
        function [dataFE, initialDataFE] = initializeMassDataFE(obj)
            dataFE = struct('Me', [], 'indexMeI', [], 'indexMeJ', []);
            initialDataFE = dataFE;
        end
        function [dataFE, initialDataFE] = initializeDataFE(obj, computePostData)
            if ~computePostData
                dataFE = struct('indexReI', [], 'Re', [], 'indexKeI', [], 'indexKeJ', [], 'Ke', []);
            else
                dataFE = struct('indexReI', [], 'Se', [], 'indexKeI', [], 'indexKeJ', [], 'Me', []);
            end
            initialDataFE = dataFE;
        end
        function dataFEMixed = initializeDataFEMixed(obj)
            dataFEMixed = struct('ReDTilde', [], 'KeDDTilde', []);
        end
        function [array, stressTensor, arrayMixed] = initializeArrayStress(obj, numberOfDofs, dimension)
            array = struct('Re', zeros(numberOfDofs, 1), 'Ke', zeros(numberOfDofs), ... % solver: computing residual and tangent
                'Se', zeros(numberOfDofs/dimension, 1), 'Me', zeros(numberOfDofs/dimension)); % postprocessing: computing stress etc.
            stressTensor = struct('FirstPK', 'Cauchy');
        end
        function dataFE = assignElementMassDataFE(obj, dataFE, Me, globalFullEdof, e)
            [globalFullEdof1, globalFullEdof2] = expandEdof(globalFullEdof(e, :));
            dataFE.indexMeI = globalFullEdof1';
            dataFE.indexMeJ = globalFullEdof2';
            dataFE.Me = Me(:);
        end
        function dataFE = assignElementDataFE(obj, dataFE, array, globalFullEdof, e, dimension, additionalFields, N, computePostData)
            if ~computePostData % residual and tangent
                dataFE.indexReI = double(globalFullEdof(e, :))';
                [globalFullEdof1, globalFullEdof2] = expandEdof(globalFullEdof(e, :));
                dataFE.indexKeI = globalFullEdof1';
                dataFE.indexKeJ = globalFullEdof2';
                dataFE.Re = array.Re;
                dataFE.Ke = array.Ke(:);
            else % compute post data (e.g. stress)
                globalFullEdof = globalFullEdof(:, 1:dimension+additionalFields:(dimension + additionalFields)*size(N, 2));
                dataFE.indexReI = double(globalFullEdof(e, :))';
                [globalFullEdof1, globalFullEdof2] = expandEdof(globalFullEdof(e, :));
                dataFE.indexKeI = globalFullEdof1';
                dataFE.indexKeJ = globalFullEdof2';
                dataFE.Se = array.Se;
                dataFE.Me = array.Me(:);
            end
        end
    end
end
