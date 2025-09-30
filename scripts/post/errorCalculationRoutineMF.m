function errorL2 = errorCalculationRoutineMF(errorL2, obj, objReference, evaluationIndex, setupObject)

%% Shape functions
N = obj.shapeFunctionObject.N_k_I;
N_k_I = obj.shapeFunctionObject.N_k_I;
dN_xi_k_I = obj.shapeFunctionObject.dN_xi_k_I;
if ~isempty(obj.mixedFEObject.shapeFunctionObject)
    M_k_I = obj.mixedFEObject.shapeFunctionObject.N_k_I;
    numberOfInternalNodes = size(M_k_I, 2);
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

% aquire material data
materialObject = obj.materialObject;
a = materialObject.a;
b = materialObject.b;
c = materialObject.c;
d = materialObject.d;

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
    saveJ = 0;
    savecReference = 0;
    saveDeltaJ = 0;
    saveLambdaC = 0;
    saveLambdaCReference = 0;
    saveDeltaLambdaC = 0;
    saveLambdaG = 0;
    saveLambdaGReference = 0;
    saveDeltaLambdaG = 0;
    saveLambdaJ = 0;
    saveLambdaJReference = 0;
    saveDeltaLambdaJ = 0;
    saveSN1 = 0;
    saveSN1Reference = 0;
    saveDeltaSN1 = 0;    
    saveSN05 = 0;
    saveSN05Reference = 0;
    saveDeltaSN05 = 0;    
end

%% element loop
numberOfElements = size(edof, 1);
for e = 1:numberOfElements
    % get element values for variables
    edN1 = qN1(edof(e, :), 1:dimension)';
    edN = obj.qN(edof(e, :), 1:dimension).';
    edN05 = 0.5*(edN + edN1);
    edN1Reference = qN1Reference(edof(e, :), 1:dimension)';
    edNReference = objReference.qN(edof(e, :), 1:dimension).';
    edN05Reference = 0.5*(edN + edN1);
        
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
            extractedJN1 = edAlphaN1(12*numberOfInternalNodes+1:13*numberOfInternalNodes).';
            extractedLambdaCN1v = edAlphaN1(13*numberOfInternalNodes+1:19*numberOfInternalNodes).';
            extractedLambdaGN1v = edAlphaN1(19*numberOfInternalNodes+1:25*numberOfInternalNodes).';
            extractedLambdaJN1 = edAlphaN1(25*numberOfInternalNodes+1:26*numberOfInternalNodes).';

            extractedCN1vReference = edAlphaN1Reference(1:6*numberOfInternalNodes).';
            extractedGN1vReference = edAlphaN1Reference(6*numberOfInternalNodes+1:12*numberOfInternalNodes).';
            extractedJN1Reference = edAlphaN1Reference(12*numberOfInternalNodes+1:13*numberOfInternalNodes).';
            extractedLambdaCN1vReference = edAlphaN1Reference(13*numberOfInternalNodes+1:19*numberOfInternalNodes).';
            extractedLambdaGN1vReference = edAlphaN1Reference(19*numberOfInternalNodes+1:25*numberOfInternalNodes).';
            extractedLambdaJN1Reference = edAlphaN1Reference(25*numberOfInternalNodes+1:26*numberOfInternalNodes).';
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
    edR = qR(edof(e, :), 1:dimension)';
    % J = edRef * dNrAll';
    % compute Jacobian
    JAll = computeJacobianForAllGausspoints(edR, dN_xi_k_I);
    %% Run through all gauss points
    for k = 1:numberOfGausspoints
        [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
        dN_X_I = computedN_X_I(dN_xi_k_I, J, k);

        Phi = edN1 * N(k, :)';
        PhiReference = edN1Reference * N(k, :)';
        deltaQ = Phi - PhiReference;
        saveQ = saveQ + Phi' * Phi * detJ * gaussWeight(k);
        saveQReference = saveQReference + PhiReference' * PhiReference * detJ * gaussWeight(k);
        saveDeltaQ = saveDeltaQ + deltaQ' * deltaQ * detJ * gaussWeight(k);
        if exist('saveDeltaC', 'var') == 1
            % C
            C = voigtToMatrix(reshape(extractedCN1v, 6, [])*M_k_I(k, :)', 'stress');
            CReference = voigtToMatrix(reshape(extractedCN1vReference, 6, [])*M_k_I(k, :)', 'stress');
            deltaC = C - CReference;
            saveC = saveC + innerProduct(C, C) * detJ * gaussWeight(k);
            saveCReference = saveCReference + innerProduct(CReference, CReference) * detJ * gaussWeight(k);
            saveDeltaC = saveDeltaC + innerProduct(deltaC, deltaC) * detJ * gaussWeight(k);
            % G
            G = voigtToMatrix(reshape(extractedGN1v, 6, [])*M_k_I(k, :)', 'stress');
            GReference = voigtToMatrix(reshape(extractedGN1vReference, 6, [])*M_k_I(k, :)', 'stress');
            deltaG = G - GReference;
            saveG = saveG + innerProduct(G, G) * detJ * gaussWeight(k);
            saveGReference = saveGReference + innerProduct(GReference, GReference) * detJ * gaussWeight(k);
            saveDeltaG = saveDeltaG + innerProduct(deltaG, deltaG) * detJ * gaussWeight(k);
            % J
            J = extractedJN1.'*M_k_I(k, :)';
            JReference = extractedJN1Reference.'*M_k_I(k, :)';
            deltaJ = J - JReference;
            saveJ = saveJ + J * J * detJ * gaussWeight(k);
            savecReference = savecReference + JReference * JReference * detJ * gaussWeight(k);
            saveDeltaJ = saveDeltaJ + deltaJ * deltaJ * detJ * gaussWeight(k);
            % LambdaC
            lambdaC = voigtToMatrix(reshape(extractedLambdaCN1v, 6, [])*M_k_I(k, :)', 'stress');
            lambdaCReference = voigtToMatrix(reshape(extractedLambdaCN1vReference, 6, [])*M_k_I(k, :)', 'stress');
            deltaLambdaC = lambdaC - lambdaCReference;
            saveLambdaC = saveLambdaC + innerProduct(lambdaC, lambdaC) * detJ * gaussWeight(k);
            saveLambdaCReference = saveLambdaCReference + innerProduct(lambdaCReference, lambdaCReference) * detJ * gaussWeight(k);
            saveDeltaLambdaC = saveDeltaLambdaC + innerProduct(deltaLambdaC, deltaLambdaC) * detJ * gaussWeight(k);
            % LambdaG
            lambdaG = voigtToMatrix(reshape(extractedLambdaGN1v, 6, [])*M_k_I(k, :)', 'stress');
            lambdaGReference = voigtToMatrix(reshape(extractedLambdaGN1vReference, 6, [])*M_k_I(k, :)', 'stress');
            deltaLambdaG = lambdaG - lambdaGReference;
            saveLambdaG = saveLambdaG + innerProduct(lambdaG, lambdaG) * detJ * gaussWeight(k);
            saveLambdaGReference = saveLambdaGReference + innerProduct(lambdaGReference, lambdaGReference) * detJ * gaussWeight(k);
            saveDeltaLambdaG = saveDeltaLambdaG + innerProduct(deltaLambdaG, deltaLambdaG) * detJ * gaussWeight(k);
            % LambdaJ
            lambdaJ = extractedLambdaJN1.'*M_k_I(k, :)';
            lambdaJReference = extractedLambdaJN1Reference.'*M_k_I(k, :)';
            deltaLambdaJ = lambdaJ - lambdaJReference;
            saveLambdaJ = saveLambdaJ + lambdaJ * lambdaJ * detJ * gaussWeight(k);
            saveLambdaJReference = saveLambdaJReference + lambdaJReference * lambdaJReference * detJ * gaussWeight(k);
            saveDeltaLambdaJ = saveDeltaLambdaJ + deltaLambdaJ * deltaLambdaJ * detJ * gaussWeight(k);
            % S
            FxN1 = edN1 * dN_X_I';
            CxN1 = FxN1.' * FxN1;
            GxN1 = 0.5 * wedge(CxN1, CxN1);
            JxN1 = det(FxN1);
            FxN1Reference = edN1Reference * dN_X_I';
            CxN1Reference = FxN1Reference.' * FxN1Reference;
            GxN1Reference = 0.5 * wedge(CxN1Reference, CxN1Reference);
            JxN1Reference = det(FxN1Reference);
            I = eye(3);
            dW_C = a * I;
            dW_G = b * I;
            dW_J = c*(J - 1) - d*J^(-1);
            dW_JReference = c*(JReference - 1) - d*JReference^(-1);
            SN1 = 2 * (dW_C + wedge(dW_G, CxN1) + 1/2*dW_J*JxN1^(-1)*GxN1);
            SN1Reference = 2 * (dW_C + wedge(dW_G, CxN1Reference) + 1/2*dW_JReference*JxN1Reference^(-1)*GxN1Reference);
            deltaSN1 = SN1 - SN1Reference;
            saveSN1 = saveSN1 + innerProduct(SN1,SN1) * detJ * gaussWeight(k);
            saveSN1Reference = saveSN1Reference + innerProduct(SN1Reference,SN1Reference) * detJ * gaussWeight(k);
            saveDeltaSN1 = saveDeltaSN1 + innerProduct(deltaSN1,deltaSN1) * detJ * gaussWeight(k);
            % 
            FxN05 = edN05 * dN_X_I';
            CxN05 = FxN05.' * FxN05;
            GxN05 = 0.5 * wedge(CxN05, CxN05);
            JxN05 = det(FxN05);            
            SN05 = 2 * (lambdaC + wedge(lambdaG, CxN05) + 1/2*lambdaJ*JxN05^(-1)*GxN05);
            FxN05Reference = edN05Reference * dN_X_I';
            CxN05Reference = FxN05Reference.' * FxN05Reference;
            GxN05Reference = 0.5 * wedge(CxN05Reference, CxN05Reference);
            JxN05Reference = det(FxN05Reference);            
            SN05Reference = 2 * (lambdaCReference + wedge(lambdaGReference, CxN05Reference) + 1/2*lambdaJReference*JxN05Reference^(-1)*GxN05Reference);     
            deltaSN05 = SN05 - SN05Reference;
            saveSN05 = saveSN05 + innerProduct(SN05,SN05) * detJ * gaussWeight(k);
            saveSN05Reference = saveSN05Reference + innerProduct(SN05Reference,SN05Reference) * detJ * gaussWeight(k);
            saveDeltaSN05 = saveDeltaSN05 + innerProduct(deltaSN05,deltaSN05) * detJ * gaussWeight(k);
        end
    end
end

%% Compute L2 norm
errorL2(evaluationIndex).phi = sqrt(saveDeltaQ/saveQReference);
if exist('saveDeltaElectroPhi', 'var') == 1
    if saveElectroPhiReference == 0
        errorL2(evaluationIndex).electroPhi = 0;
    else
        errorL2(evaluationIndex).electroPhi = sqrt(saveDeltaElectroPhi/saveElectroPhiReference);
    end
end
if exist('saveDeltaC', 'var') == 1
    errorL2(evaluationIndex).C = sqrt(saveDeltaC/saveCReference);
    errorL2(evaluationIndex).G = sqrt(saveDeltaG/saveGReference);
    errorL2(evaluationIndex).J = sqrt(saveDeltaJ/savecReference);
    errorL2(evaluationIndex).lambdaC = sqrt(saveDeltaLambdaC/saveLambdaCReference);
    errorL2(evaluationIndex).lambdaG = sqrt(saveDeltaLambdaG/saveLambdaGReference);
    errorL2(evaluationIndex).lambda = sqrt(saveDeltaLambdaJ/saveLambdaJReference);
    errorL2(evaluationIndex).SN05 = sqrt(saveDeltaSN05/saveSN05Reference);    
    errorL2(evaluationIndex).SN1 = sqrt(saveDeltaSN1/saveSN1Reference);    
end
end
