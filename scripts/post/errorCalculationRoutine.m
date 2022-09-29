function errorL2 = errorCalculationRoutine(errorL2, obj, objReference, evaluationIndex)

%% Shape functions
N = obj.shapeFunctionObject.N;
dNrAll = obj.shapeFunctionObject.dNr;
if ~isempty(obj.mixedFEObject.shapeFunctionObject)
    M = obj.mixedFEObject.shapeFunctionObject.N;
    numberOfInternalNodes = size(M, 2);
end
gaussWeight = obj.shapeFunctionObject.gaussWeight;
numberOfGausspoints = obj.shapeFunctionObject.numberOfGausspoints;
dimension = obj.dimension;
edof = obj.meshObject.edof;

% variables
qR = obj.qR;
qN1 = obj.qN1;
qN1Reference = objReference.qN1;
alphaN1 = obj.mixedFEObject.qN1;
alphaN1Reference = objReference.mixedFEObject.qN1;

% initialize save variables
saveQ = 0;
saveQReference = 0;
saveDeltaQ = 0;
if isa(obj, 'solidElectroClass') || isa(obj, 'solidElectroThermoClass')
    saveElectroPhi = 0;
    saveElectroPhiReference = 0;
    saveDeltaElectroPhi = 0;

    saveD = 0;
    saveDReference = 0;
    saveDeltaD = 0;
end
if isa(obj, 'solidThermoClass') || isa(obj, 'solidElectroThermoClass')
    saveTheta = 0;
    saveThetaReference = 0;
    saveDeltaTheta = 0;
end
if strcmp(obj.elementDisplacementType, 'mixedSC')
    saveC = 0;
    saveCReference = 0;
    saveDeltaC = 0;
    saveG = 0;
    saveGReference = 0;
    saveDeltaG = 0;
    savec = 0;
    savecReference = 0;
    saveDeltac = 0;
    saveLambdaC = 0;
    saveLambdaCReference = 0;
    saveDeltaLambdaC = 0;
    saveLambdaG = 0;
    saveLambdaGReference = 0;
    saveDeltaLambdaG = 0;
    saveLambdac = 0;
    saveLambdacReference = 0;
    saveDeltaLambdac = 0;
end

%% element loop
numberOfElements = size(edof, 1);
for e = 1:numberOfElements
    % get element values for variables
    edN1 = qN1(edof(e, :), 1:dimension)';
    edN1Reference = qN1Reference(edof(e, :), 1:dimension)';

    if isa(obj, 'solidThermoClass')
        thetaN1 = qN1(edof(e, :), dimension+1)';
        thetaN1Reference = qN1Reference(edof(e, :), dimension+1)';
    elseif isa(obj, 'solidElectroClass') || isa(obj, 'solidElectroThermoClass')
        phiN1 = qN1(edof(e, :), dimension+1)';
        phiN1Reference = qN1Reference(edof(e, :), dimension+1)';

        if isa(obj, 'solidElectroThermoClass')
            thetaN1 = qN1(edof(e, :), dimension+2)';
            thetaN1Reference = qN1Reference(edof(e, :), dimension+2)';
        end
    end
    if strcmp(obj.elementDisplacementType, 'mixedSC')
        edAlphaN1 = alphaN1(e, :);
        edAlphaN1Reference = alphaN1Reference(e, :);

        if isa(obj, 'solidClass') || isa(obj, 'solidThermoClass')
            extractedCN1v = edAlphaN1(1:6*numberOfInternalNodes).';
            extractedGN1v = edAlphaN1(6*numberOfInternalNodes+1:12*numberOfInternalNodes).';
            extractedcN1 = edAlphaN1(12*numberOfInternalNodes+1:13*numberOfInternalNodes).';
            extractedLambdaCN1v = edAlphaN1(13*numberOfInternalNodes+1:19*numberOfInternalNodes).';
            extractedLambdaGN1v = edAlphaN1(19*numberOfInternalNodes+1:25*numberOfInternalNodes).';
            extractedLambdacN1 = edAlphaN1(25*numberOfInternalNodes+1:26*numberOfInternalNodes).';

            extractedCN1vReference = edAlphaN1Reference(1:6*numberOfInternalNodes).';
            extractedGN1vReference = edAlphaN1Reference(6*numberOfInternalNodes+1:12*numberOfInternalNodes).';
            extractedcN1Reference = edAlphaN1Reference(12*numberOfInternalNodes+1:13*numberOfInternalNodes).';
            extractedLambdaCN1vReference = edAlphaN1Reference(13*numberOfInternalNodes+1:19*numberOfInternalNodes).';
            extractedLambdaGN1vReference = edAlphaN1Reference(19*numberOfInternalNodes+1:25*numberOfInternalNodes).';
            extractedLambdacN1Reference = edAlphaN1Reference(25*numberOfInternalNodes+1:26*numberOfInternalNodes).';
        elseif isa(obj, 'solidElectroClass') || isa(obj, 'solidElectroThermoClass')
            extractedDN1 = edAlphaN1(1:3*numberOfInternalNodes).';
            extractedCN1v = edAlphaN1(3*numberOfInternalNodes+1:9*numberOfInternalNodes).';
            extractedGN1v = edAlphaN1(9*numberOfInternalNodes+1:15*numberOfInternalNodes).';
            extractedcN1 = edAlphaN1(15*numberOfInternalNodes+1:16*numberOfInternalNodes).';
            extractedLambdaCN1v = edAlphaN1(16*numberOfInternalNodes+1:22*numberOfInternalNodes).';
            extractedLambdaGN1v = edAlphaN1(22*numberOfInternalNodes+1:28*numberOfInternalNodes).';
            extractedLambdacN1 = edAlphaN1(28*numberOfInternalNodes+1:29*numberOfInternalNodes).';

            extractedDN1Reference = edAlphaN1Reference(1:3*numberOfInternalNodes).';
            extractedCN1vReference = edAlphaN1Reference(3*numberOfInternalNodes+1:9*numberOfInternalNodes).';
            extractedGN1vReference = edAlphaN1Reference(9*numberOfInternalNodes+1:15*numberOfInternalNodes).';
            extractedcN1Reference = edAlphaN1Reference(15*numberOfInternalNodes+1:16*numberOfInternalNodes).';
            extractedLambdaCN1vReference = edAlphaN1Reference(16*numberOfInternalNodes+1:22*numberOfInternalNodes).';
            extractedLambdaGN1vReference = edAlphaN1Reference(22*numberOfInternalNodes+1:28*numberOfInternalNodes).';
            extractedLambdacN1Reference = edAlphaN1Reference(28*numberOfInternalNodes+1:29*numberOfInternalNodes).';
        end
    end

    % Jacobi-Matrix
    edRef = qR(edof(e, :), 1:dimension)';
    J = edRef * dNrAll';

    %% Run through all gauss points
    for k = 1:numberOfGausspoints
        indx = dimension * k - (dimension - 1):dimension * k;
        detJ = det(J(:, indx)');
        if detJ < 10 * eps
            error('Jacobi determinant equal or less than zero.')
        end

        Phi = edN1 * N(k, :)';
        PhiReference = edN1Reference * N(k, :)';
        deltaQ = Phi - PhiReference;
        saveQ = saveQ + Phi' * Phi * detJ * gaussWeight(k);
        saveQReference = saveQReference + PhiReference' * PhiReference * detJ * gaussWeight(k);
        saveDeltaQ = saveDeltaQ + deltaQ' * deltaQ * detJ * gaussWeight(k);
        
        
        if exist('saveDeltaElectroPhi', 'var') == 1
            % electro phi
            electroPhi = phiN1 * N(k, :)';
            electroPhiReference = phiN1Reference * N(k, :)';
            deltaElectroPhi = electroPhi - electroPhiReference;
            saveElectroPhi = saveElectroPhi + electroPhi * electroPhi * detJ * gaussWeight(k);
            saveElectroPhiReference = saveElectroPhiReference + electroPhiReference * electroPhiReference * detJ * gaussWeight(k);
            saveDeltaElectroPhi = saveDeltaElectroPhi + deltaElectroPhi * deltaElectroPhi * detJ * gaussWeight(k);
        end
        if exist('saveDeltaD', 'var') == 1
            % D
            D = reshape(extractedDN1, 3, []) * M(k, :)';
            DReference = reshape(extractedDN1Reference, 3, []) * M(k, :)';
            deltaD = D - DReference;
            saveD = saveD + D' * D * detJ * gaussWeight(k);
            saveDReference = saveDReference + DReference' * DReference * detJ * gaussWeight(k);
            saveDeltaD = saveDeltaD + deltaD' * deltaD * detJ * gaussWeight(k);
        end
        if exist('saveDeltaTheta', 'var') == 1
            % theta
            theta = thetaN1 * N(k, :)';
            thetaReference = thetaN1Reference * N(k, :)';
            deltaTheta = theta - thetaReference;
            saveTheta = saveTheta + theta * theta * detJ * gaussWeight(k);
            saveThetaReference = saveThetaReference + thetaReference * thetaReference * detJ * gaussWeight(k);
            saveDeltaTheta = saveDeltaTheta + deltaTheta * deltaTheta * detJ * gaussWeight(k);
        end
        if exist('saveDeltaC', 'var') == 1
            % C
            C = voigtToMatrix(reshape(extractedCN1v, 6, [])*M(k, :)', 'stress');
            CReference = voigtToMatrix(reshape(extractedCN1vReference, 6, [])*M(k, :)', 'stress');
            deltaC = C - CReference;
            saveC = saveC + innerProduct(C, C) * detJ * gaussWeight(k);
            saveCReference = saveCReference + innerProduct(CReference, CReference) * detJ * gaussWeight(k);
            saveDeltaC = saveDeltaC + innerProduct(deltaC, deltaC) * detJ * gaussWeight(k);
            % G
            G = voigtToMatrix(reshape(extractedGN1v, 6, [])*M(k, :)', 'stress');
            GReference = voigtToMatrix(reshape(extractedGN1vReference, 6, [])*M(k, :)', 'stress');
            deltaG = G - GReference;
            saveG = saveG + innerProduct(G, G) * detJ * gaussWeight(k);
            saveGReference = saveGReference + innerProduct(GReference, GReference) * detJ * gaussWeight(k);
            saveDeltaG = saveDeltaG + innerProduct(deltaG, deltaG) * detJ * gaussWeight(k);
            % c
            c = extractedcN1.'*M(k, :)';
            cReference = extractedcN1Reference.'*M(k, :)';
            deltac = c - cReference;
            savec = savec + c * c * detJ * gaussWeight(k);
            savecReference = savecReference + cReference * cReference * detJ * gaussWeight(k);
            saveDeltac = saveDeltac + deltac * deltac * detJ * gaussWeight(k);
            % LambdaC
            lambdaC = voigtToMatrix(reshape(extractedLambdaCN1v, 6, [])*M(k, :)', 'stress');
            lambdaCReference = voigtToMatrix(reshape(extractedLambdaCN1vReference, 6, [])*M(k, :)', 'stress');
            deltaLambdaC = lambdaC - lambdaCReference;
            saveLambdaC = saveLambdaC + innerProduct(lambdaC, lambdaC) * detJ * gaussWeight(k);
            saveLambdaCReference = saveLambdaCReference + innerProduct(lambdaCReference, lambdaCReference) * detJ * gaussWeight(k);
            saveDeltaLambdaC = saveDeltaLambdaC + innerProduct(deltaLambdaC, deltaLambdaC) * detJ * gaussWeight(k);
            % LambdaG
            lambdaG = voigtToMatrix(reshape(extractedLambdaGN1v, 6, [])*M(k, :)', 'stress');
            lambdaGReference = voigtToMatrix(reshape(extractedLambdaGN1vReference, 6, [])*M(k, :)', 'stress');
            deltaLambdaG = lambdaG - lambdaGReference;
            saveLambdaG = saveLambdaG + innerProduct(lambdaG, lambdaG) * detJ * gaussWeight(k);
            saveLambdaGReference = saveLambdaGReference + innerProduct(lambdaGReference, lambdaGReference) * detJ * gaussWeight(k);
            saveDeltaLambdaG = saveDeltaLambdaG + innerProduct(deltaLambdaG, deltaLambdaG) * detJ * gaussWeight(k);
            % Lambdac
            lambdac = extractedLambdacN1.'*M(k, :)';
            lambdacReference = extractedLambdacN1Reference.'*M(k, :)';
            deltaLambdac = lambdac - lambdacReference;
            saveLambdac = saveLambdac + lambdac * lambdac * detJ * gaussWeight(k);
            saveLambdacReference = saveLambdacReference + lambdacReference * lambdacReference * detJ * gaussWeight(k);
            saveDeltaLambdac = saveDeltaLambdac + deltaLambdac * deltaLambdac * detJ * gaussWeight(k);
        end

    end
end

%% Compute L2 norm
errorL2(evaluationIndex).Q = sqrt(saveDeltaQ/saveQReference);
if exist('saveDeltaElectroPhi', 'var') == 1
    if saveElectroPhiReference == 0
        errorL2(evaluationIndex).electroPhi = 0;
    else
        errorL2(evaluationIndex).electroPhi = sqrt(saveDeltaElectroPhi/saveElectroPhiReference);
    end
end
if exist('saveDeltaD', 'var') == 1
    if saveDReference == 0
        errorL2(evaluationIndex).D = 0;
    else
        errorL2(evaluationIndex).D = sqrt(saveDeltaD/saveDReference);
    end
end
if exist('saveDeltaTheta', 'var') == 1
    errorL2(evaluationIndex).theta = sqrt(saveDeltaTheta/saveThetaReference);
end
if exist('saveDeltaC', 'var') == 1
    errorL2(evaluationIndex).C = sqrt(saveDeltaC/saveCReference);
    errorL2(evaluationIndex).G = sqrt(saveDeltaG/saveGReference);
    errorL2(evaluationIndex).c = sqrt(saveDeltac/savecReference);
    errorL2(evaluationIndex).lambdaC = sqrt(saveDeltaLambdaC/saveLambdaCReference);
    errorL2(evaluationIndex).lambdaG = sqrt(saveDeltaLambdaG/saveLambdaGReference);
    errorL2(evaluationIndex).lambdac = sqrt(saveDeltaLambdac/saveLambdacReference);
end
end
