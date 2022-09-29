classdef saveClass < matlab.mixin.Copyable
    properties
        saveData = false;
        fileName = 'simulationData';

        reduceFileSize = true;
        saveForTimeSteps = [];
        saveDataInOneFile = true;

        exportLogFile = false;
    end
    methods
        function initialize(obj, setupObject)
            % begin writing to logfile
            generatedFilesFolderInfo = what('generatedFiles');
            if obj.exportLogFile
                filenameLogFile = strcat(generatedFilesFolderInfo.path, filesep, 'logFiles', filesep, obj.fileName, 'Log.txt');
                if isfile(filenameLogFile)
                    delete(filenameLogFile);
                end
                diary(filenameLogFile);
            end

            % initialize saveObject
            if obj.saveData
                if obj.saveDataInOneFile
                    % delete .mat file if it already exists
                    filename = strcat(generatedFilesFolderInfo.path, filesep, 'matFiles', filesep, obj.fileName, '.mat');
                    if isfile(filename)
                        delete(filename);
                    end
                end

                % set default configuration (save data for all time steps)
                if isempty(obj.saveForTimeSteps)
                    obj.saveForTimeSteps = 1:setupObject.totalTimeSteps;
                end
            end
        end
        function saveDataAsMatFile(obj, setupObject, dofObject, timeStep)
            if obj.saveData
                if ismember(timeStep, obj.saveForTimeSteps)
                    % reduce file size by deleting irrelevant entries from dofObject before saving
                    if obj.reduceFileSize
                        listContinuumObjects = dofObject.listContinuumObjects;
                        for jj = 1:size(listContinuumObjects, 2)
                            if isa(listContinuumObjects{jj}, 'solidSuperClass')
                                continuumObject = listContinuumObjects{jj};
                                dataFE = continuumObject.storageFEObject.dataFE;
                                mixedDataFE = continuumObject.mixedFEObject.dataFE;
                                continuumObject.storageFEObject.dataFE = struct();
                                continuumObject.mixedFEObject.dataFE = struct();
                            end
                        end
                    end

                    % save simulation data in .mat file(s)
                    generatedFilesFolderInfo = what('generatedFiles');
                    filename = strcat(generatedFilesFolderInfo.path, filesep, 'matFiles', filesep, obj.fileName);
                    if obj.saveDataInOneFile
                        filename = strcat(filename, '.mat');
                        eval(strcat('setupObject', num2str(timeStep), ' = setupObject;'));
                        eval(strcat('dofObject', num2str(timeStep), ' = dofObject;'));
                        if timeStep == min(obj.saveForTimeSteps)
                            % create new .mat file and save simulation data
                            save(filename, strcat('setupObject', num2str(timeStep)), strcat('dofObject', num2str(timeStep)));
                        else
                            % append simulation data to existing .mat file
                            save(filename, strcat('setupObject', num2str(timeStep)), strcat('dofObject', num2str(timeStep)), '-append');
                        end
                    else
                        filename = strcat(filename, num2str(timeStep), '.mat');
                        save(filename, 'setupObject', 'dofObject');
                    end

                    % another implementation idea using matfile (slower)
                    % m = matfile(strcat('matFiles/', setupObject.saveObject.fileName),'Writable',true);
                    % m.saveJournal(timeStep,1) = struct('setupObject', setupObject, 'dofObject', dofObject);

                    % reset dofObject if entries were deleted
                    if obj.reduceFileSize
                        for jj = 1:size(listContinuumObjects, 2)
                            if isa(listContinuumObjects{jj}, 'solidSuperClass')
                                continuumObject = listContinuumObjects{jj};
                                continuumObject.storageFEObject.dataFE = dataFE;
                                continuumObject.mixedFEObject.dataFE = mixedDataFE;
                            end
                        end
                    end
                end
            end
        end

        function terminate(obj)
            % end writing to logfile
            if obj.exportLogFile
                diary off;
            end
        end

        %% set and get methods
        function set.saveData(obj, val)
            assert(islogical(val), 'saveData must be logical')
            obj.saveData = val;
        end

        function set.fileName(obj, val)
            assert(ischar(val), 'fileName must be of type string!');
            obj.fileName = val;
        end

        function set.reduceFileSize(obj, val)
            assert(islogical(val), 'reduceFileSize must be logical')
            obj.reduceFileSize = val;
        end

        function set.saveForTimeSteps(obj, val)
            assert(isvector(val), 'saveForTimeSteps must be a vector')
            obj.saveForTimeSteps = val;
        end

        function set.saveDataInOneFile(obj, val)
            assert(islogical(val), 'saveDataInOneFile must be logical')
            obj.saveDataInOneFile = val;
        end

        function set.exportLogFile(obj, val)
            assert(islogical(val), 'exportLogFile must be logical')
            obj.exportLogFile = val;
        end
    end
end
