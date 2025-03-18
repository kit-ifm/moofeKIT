function dataFEContinuumObject = callElements(continuumObject, setupObject, requiredData)
%CALLELEMENTS Loop over all finite elements
%
%   CALL
%   dataFEContinuumObject = callElements(continuumObject, setupObject, requiredData)
%   continuumObject:    continuumObject (e.g. object of solidClass, solidThermoClass, ...)
%   setupObject:        setupObject
%   requiredData:       string that specifies the required data ('residualAndTangent', 'massMatrix', 'postData')

%% FIXME: unify (for all classes)
computePostData = false;
if strcmp(requiredData, 'postData')
    computePostData = true;
end

% assemble elementName for function call
if isa(continuumObject, 'solidSuperClass')
    elementDisplacementType = continuumObject.elementDisplacementType;
    elementMaterialName = continuumObject.materialObject.name;
    elementNameAdditionalSpecification = continuumObject.elementNameAdditionalSpecification;
    if continuumObject.dimension == 2
        switch continuumObject.elementDisplacementType
            case {'incompatibleModesWilson', 'incompatibleModesTaylor'}
                elementDisplacementType = 'incompatibleModes';
        end
        switch continuumObject.materialObject.name
            case {'HookeESZ', 'HookeEVZ'}
                elementMaterialName = 'Hooke2D';
            case {'SaintVenantESZ', 'SaintVenantEVZ'}
                elementMaterialName = 'SaintVenant2D';
            case {'NeoHookeESZ', 'NeoHookeEVZ'}
                elementMaterialName = 'NeoHooke2D';
        end
    end
    if isa(continuumObject, 'beamClass')
        elementNameAdditionalSpecification = strcat(elementNameAdditionalSpecification, continuumObject.theory);
    end
    elementName = strcat(elementDisplacementType, elementNameAdditionalSpecification, elementMaterialName, setupObject.integrator);

elseif isa(continuumObject, 'neumannClass')
    elementName = continuumObject.loadType;
elseif isa(continuumObject, 'dirichletLagrangeClass')
    elementName = strcat('nodal', setupObject.integrator);
elseif isa(continuumObject, 'bodyForceClass')
    elementName = strcat(continuumObject.typeOfLoad, setupObject.integrator);
else
    error('element name not available!')
end

if isa(continuumObject, 'solidSuperClass') || isa(continuumObject, 'neumannClass') || isa(continuumObject, 'bodyForceClass')

    %% call single element
    % load objects
    shapeFunctionObject = continuumObject.shapeFunctionObject;
    meshObject = continuumObject.meshObject;
    storageFEObject = continuumObject.storageFEObject;

    % number of degrees of freedom and more
    N_k_I = shapeFunctionObject.N_k_I;
    edof = meshObject.edof;
    numberOfNodes = shapeFunctionObject.numberOfNodes;
    if isa(continuumObject, 'solidSuperClass')
        dimension = continuumObject.dimension;
        mixedFEObject = continuumObject.mixedFEObject;
        numericalTangentObject = continuumObject.numericalTangentObject;

        globalEdofContinuumObject = meshObject.globalFullEdof;
        if isa(continuumObject, 'plateClass') || isa(continuumObject, 'stringClass') || isa(continuumObject, 'beamClass')
            numberOfPrimaryDofsPerNode = continuumObject.numberOfDofsPerNode;
        else
            numberOfPrimaryDofsPerNode = dimension;
        end
        
        numberOfAdditionalNonMixedFields = continuumObject.additionalFields;
        numberOfMixedFields = mixedFEObject.numberOfFields;
        dofsPerMixedField = mixedFEObject.dofsPerField;
        if ~isempty(dofsPerMixedField)
            temp = vertcat(dofsPerMixedField{:});
            dofsPerMixedField = temp(:,1);
        else
            dofsPerMixedField = [];
        end

        if ~isempty(continuumObject.dofsPerAdditionalField)
            numberOfAdditionalDofs = continuumObject.dofsPerAdditionalField * numberOfNodes;
        else
            numberOfAdditionalDofs = [];
        end

    elseif isa(continuumObject, 'neumannClass')
        globalEdofContinuumObject = continuumObject.globalEdof;
        numberOfAdditionalDofsPerNode = sum(continuumObject.masterObject.dofsPerAdditionalField);
        numberOfPrimaryDofsPerNode = size(continuumObject.masterObject.qR, 2) - numberOfAdditionalDofsPerNode;
        numberOfAdditionalNonMixedFields = 0;
        numberOfMixedFields = 0;
        dofsPerMixedField = 0;
        dimension = 0;
        numberOfDofsNeumann = size(globalEdofContinuumObject, 2);
    elseif isa(continuumObject, 'bodyForceClass')
        dimension = continuumObject.dimension;
        globalEdofContinuumObject = continuumObject.meshObject.globalFullEdof;
        numberOfPrimaryDofsPerNode = 0;
        numberOfAdditionalNonMixedFields = 0;
        numberOfMixedFields = 0;
        dofsPerMixedField = 0;
        numberOfDofsBodyForce = dimension * size(meshObject.edof, 2);        
    end
    if length(numberOfNodes) == 1
        numberOfPrimaryDofs = numberOfPrimaryDofsPerNode * numberOfNodes;
    else
        numberOfPrimaryDofs =  sum(numberOfNodes(1:numberOfPrimaryDofsPerNode));
    end
    numberOfNonMixedFields = 1 + numberOfAdditionalNonMixedFields;
    totalNumberOfFields = numberOfNonMixedFields + numberOfMixedFields;

    % save dofs per field
    if isa(continuumObject, 'solidSuperClass')
        dofsData = [numberOfPrimaryDofs; numberOfAdditionalDofs; dofsPerMixedField];
    elseif isa(continuumObject, 'neumannClass')
        dofsData = numberOfDofsNeumann;
    elseif isa(continuumObject, 'bodyForceClass')
        dofsData = numberOfDofsBodyForce;
    end

    % set initial residual and tangent values
    rDataInitial = cell(totalNumberOfFields, 1);
    kDataInitial = cell(totalNumberOfFields, totalNumberOfFields);
    for ii = 1:totalNumberOfFields
        rDataInitial{ii, 1} = zeros(dofsData(ii), 1);
        for jj = 1:totalNumberOfFields
            kDataInitial{ii, jj} = zeros(dofsData(ii), dofsData(jj));
        end
    end

    % numericalTangent information
    numericalTangentInfo = struct('computeNumericalTangent', false, 'showDifferences', false, 'type', '', 'epsilon', 0);
    if isa(continuumObject, 'solidSuperClass')
        numericalTangentInfo.computeNumericalTangent = numericalTangentObject.computeNumericalTangent;
        numericalTangentInfo.showDifferences = numericalTangentObject.showDifferences;
        numericalTangentInfo.type = numericalTangentObject.type;
        numericalTangentInfo.epsilon = numericalTangentObject.epsilon;
    end

    % element loop
    [dataFEContinuumObject, initialDataFE] = storageFEObject.initializeDataFE(requiredData);
    if ~strcmp(requiredData, 'massMatrix')
        dataFEMixed = storageFEObject.initializeDataFEMixed;
    end
    numberOfElements = size(globalEdofContinuumObject, 1);
    % parfor e = 1:numberOfElements
    for e = 1:numberOfElements
        dofs = struct();
        if isa(continuumObject, 'solidSuperClass')
            dofs.edN1 = continuumObject.qN1(edof(e, :), 1:numberOfPrimaryDofsPerNode)';
            if isa(continuumObject, 'beamClass')
                dofs.phiN1 = continuumObject.qN1(edof(e, :), numberOfPrimaryDofsPerNode+1)'; %dimension+1
            elseif isa(continuumObject, 'solidElectroClass')
                dofs.phiN1 = continuumObject.qN1(edof(e, :), dimension+1)';
            elseif isa(continuumObject, 'solidThermoClass')
                dofs.thetaN1 = continuumObject.qN1(edof(e, :), dimension+1)';
            elseif isa(continuumObject, 'solidVelocityClass')
                dofs.vN1 = continuumObject.qN1(edof(e, :), dimension+1:end)';
            elseif isa(continuumObject, 'solidElectroThermoClass')
                dofs.phiN1 = continuumObject.qN1(edof(e, :), dimension+1)';
                dofs.thetaN1 = continuumObject.qN1(edof(e, :), dimension+2)';
            end

            edAlphaN1 = continuumObject.mixedFEObject.qN1;
            if size(edAlphaN1, 1) > 0
                dofs.edAlphaN1 = edAlphaN1(e, :);
            end
        end

        if computePostData
            [array, stressTensor] = storageFEObject.initializeArrayStress(numberOfPrimaryDofs, numberOfNodes);
        else
            array = struct();
            stressTensor = struct();
        end
        if strcmp(requiredData, 'massMatrix')
            if isa(continuumObject,'beamClass')
                array = beamMassMatrixElement(continuumObject, setupObject, e);
            else
                array = massMatrixElement(continuumObject, setupObject, e);
            end
            
        else
            [rData, kData, elementEnergy, array] = feval(elementName, continuumObject, setupObject, computePostData, e, rDataInitial, kDataInitial, dofs, array, stressTensor, false);
        end

        if strcmp(requiredData, 'residualAndTangent')
            % numerical tangent
            if numericalTangentInfo.computeNumericalTangent
                namesOfFields = fieldnames(dofs);
                numberOfFields = length(namesOfFields);
                fieldNameCell = cell(totalNumberOfFields, 1);
                fieldNameCell(1:numberOfFields) = namesOfFields;

                numericalTangentDifferences = cell(totalNumberOfFields);
                numericalTangentDifferencesNames = cell(totalNumberOfFields, 1);

                totalNumberOfDofs = sum(dofsData);
                numberOfNonMixedDofs = sum(dofsData(1:numberOfNonMixedFields));
                kDataNTMatrix = zeros(totalNumberOfDofs);
                %             parfor ii=1:totalNumberOfDofs
                for ii = 1:totalNumberOfDofs
                    if ii <= numberOfNonMixedDofs
                        indexFieldNameCell = 1;
                        for jj = 1:numberOfFields
                            if ii <= sum(dofsData(1:jj))
                                indexFieldNameCell = jj;
                                break;
                            end
                        end
                        indexDof = ii - sum(dofsData(1:indexFieldNameCell-1));
                    end
                    if any(strcmpi(numericalTangentInfo.type, {'standard', 'centralDifferences', 'complex'}))
                        dofsNTxx = dofs;
                        if ii <= numberOfNonMixedDofs
                            dofsNTxx.(fieldNameCell{indexFieldNameCell})(indexDof) = dofsNTxx.(fieldNameCell{indexFieldNameCell})(indexDof) + numericalTangentInfo.epsilon;
                        else
                            dofsNTxx.edAlphaN1(ii-numberOfNonMixedDofs) = dofsNTxx.edAlphaN1(ii-numberOfNonMixedDofs) + numericalTangentInfo.epsilon;
                        end
                        [rxxDataForward, ~] = feval(elementName, continuumObject, setupObject, computePostData, e, rDataInitial, kDataInitial, dofsNTxx, struct(), struct(), true);
                    end
                    if strcmpi(numericalTangentInfo.type, 'centralDifferences')
                        dofsNTxx = dofs;
                        if ii <= numberOfNonMixedDofs
                            dofsNTxx.(fieldNameCell{indexFieldNameCell})(indexDof) = dofsNTxx.(fieldNameCell{indexFieldNameCell})(indexDof) - numericalTangentInfo.epsilon;
                        else
                            dofsNTxx.edAlphaN1(ii-numberOfNonMixedDofs) = dofsNTxx.edAlphaN1(ii-numberOfNonMixedDofs) - numericalTangentInfo.epsilon;
                        end
                        [rxxDataBackward, ~] = feval(elementName, continuumObject, setupObject, computePostData, e, rDataInitial, kDataInitial, dofsNTxx, struct(), struct(), true);
                    end

                    if strcmpi(numericalTangentInfo.type, 'standard')
                        kDataNTMatrix(:, ii) = (cell2mat(rxxDataForward) - cell2mat(rData)) / numericalTangentInfo.epsilon;
                    elseif strcmpi(numericalTangentInfo.type, 'centralDifferences')
                        kDataNTMatrix(:, ii) = (cell2mat(rxxDataForward) - cell2mat(rxxDataBackward)) / (2 * numericalTangentInfo.epsilon);
                    elseif strcmpi(numericalTangentInfo.type, 'complex')
                        kDataNTMatrix(:, ii) = imag(cell2mat(rxxDataForward)) / imag(numericalTangentInfo.epsilon);
                    end
                end
                kDataNT = mat2cell(kDataNTMatrix, dofsData, dofsData);
                if numericalTangentInfo.showDifferences
                    for ii = 1:totalNumberOfFields
                        for jj = 1:totalNumberOfFields
                            numericalTangentDifferences{ii, jj} = max(max(abs(kData{ii, jj}-kDataNT{ii, jj})));
                            numericalTangentDifferencesNames{ii} = strcat(num2str(ii), '.');
                        end
                    end
                    numericalTangentDifferencesTable = cell2table(numericalTangentDifferences, 'VariableNames', numericalTangentDifferencesNames, 'RowNames', numericalTangentDifferencesNames);
                    disp(numericalTangentDifferencesTable);
                end
                kData = kDataNT;
            end

            % compute residuals and tangents / static condensation
            RD = cell2mat(rData(1:numberOfNonMixedFields));
            RA = cell2mat(rData(numberOfNonMixedFields+1:totalNumberOfFields));
            KDD = cell2mat(kData(1:numberOfNonMixedFields, 1:numberOfNonMixedFields));
            KDA = cell2mat(kData(1:numberOfNonMixedFields, numberOfNonMixedFields+1:totalNumberOfFields));
            KAD = cell2mat(kData(numberOfNonMixedFields+1:totalNumberOfFields, 1:numberOfNonMixedFields));
            KAA = cell2mat(kData(numberOfNonMixedFields+1:totalNumberOfFields, numberOfNonMixedFields+1:totalNumberOfFields));

            if isa(continuumObject, 'solidSuperClass')
                if isprop(continuumObject, 'permutationMatrix')
                    P = continuumObject.permutationMatrix;
                    [array, elementArrayMixed] = mixedFEFunction(continuumObject.mixedFEObject, setupObject, RD, RA, KDD, KDA, KAD, KAA, P);
                    dataFEMixed(e) = elementArrayMixed;
                else
                    [array, elementArrayMixed] = mixedFEFunction(continuumObject.mixedFEObject, setupObject, RD, RA, KDD, KDA, KAD, KAA);
                    dataFEMixed(e) = elementArrayMixed;
                end
            elseif isa(continuumObject, 'neumannClass') || isa(continuumObject, 'bodyForceClass')
                array.Re = RD;
                array.Ke = KDD;
            end

            % store energy
            globalEnergy(e) = elementEnergy;
        end

        % store FE Data
        elementDataFE = storageFEObject.assignElementDataFE(initialDataFE, array, globalEdofContinuumObject, e, numberOfPrimaryDofsPerNode, numberOfAdditionalNonMixedFields, N_k_I, requiredData);
        dataFEContinuumObject(e) = elementDataFE;
    end

    if strcmp(requiredData, 'residualAndTangent')
        if isa(continuumObject, 'solidSuperClass')
            mixedFEObject.dataFE = dataFEMixed;
        end
        
        % save energy data for each time step
        fieldnamesEnergy = fieldnames(globalEnergy);
        if setupObject.timeStep == 1 && setupObject.newton.step(setupObject.timeStep) == 1
            for ii = 1:size(fieldnamesEnergy, 1)
                continuumObject.ePot(1).(fieldnamesEnergy{ii}) = sum([globalEnergy.(fieldnamesEnergy{ii})]);
            end
        end
        for ii = 1:size(fieldnamesEnergy, 1)
            continuumObject.ePot(setupObject.timeStep+1).(fieldnamesEnergy{ii}) = sum([globalEnergy.(fieldnamesEnergy{ii})]);
        end
        clear globalEnergy;
    end

    if isa(continuumObject, 'neumannClass') || isa(continuumObject, 'bodyForceClass')
        if ~computePostData
            continuumObject.storageFEObject.dataFE = dataFEContinuumObject;
        end
    end
else
    feval(elementName, continuumObject, setupObject, computePostData)
end
end
