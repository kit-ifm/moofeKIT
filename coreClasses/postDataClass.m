classdef postDataClass < matlab.mixin.Copyable
    properties
        energyJournal = struct('time', [], 'EKin', [], 'EKinDifference', []);
        storeStateFlag = false;
        stateJournal = struct('time', [], 'position', [], 'velocity', []);
    end

    methods
        function obj = postDataClass
            % constructor
        end
        function initialize(obj, dofObject, M)
            vN = dofObject.vN;
            updateEnergyJournal(obj, 1, 0, vN, vN, M);
            updateMomenta(obj, dofObject, 1, vN, M);
        end
        function update(obj, dofObject, timeStep, time, M)
            vN = dofObject.vN;
            vN1 = dofObject.vN1;
            indexStructure = timeStep + 1;
            updateEnergyJournal(obj, indexStructure, time, vN, vN1, M);
            updateMomenta(obj, dofObject, indexStructure, vN1, M);
            if obj.storeStateFlag
                updateStateJournal(obj, dofObject, indexStructure, time, vN, vN1);
            end
        end
        function updateEnergyJournal(obj, index, time, vN, vN1, M)
            obj.energyJournal(index).time = time;
            EKinN = 1 / 2 * vN' * M * vN;
            EKinN1 = 1 / 2 * vN1' * M * vN1;
            obj.energyJournal(index).EKin = EKinN1;
            obj.energyJournal(index).EKinDifference = EKinN1 - EKinN;
        end
        function updateMomenta(obj, dofObject, index, v, M)
            dofObject.L = M * v;
            for ii = 1:dofObject.numberOfContinuumObjects
                continuumObject = dofObject.listContinuumObjects{ii};
                if isprop(continuumObject, 'momenta')
                    if isa(continuumObject,"beamClass")
                        continuumObject.L = dofObject.L(continuumObject.meshObject.globalNodesDof);
                    else
                        continuumObject.L = dofObject.L(continuumObject.meshObject.globalNodesDof(:, 1:continuumObject.dimension));
                        continuumObject.momenta(index).L = sum(continuumObject.L, 1);
                    end
                    continuumObject.momenta(index).J = continuumObject.J;
                end
            end
            0;
        end
        function updateStateJournal(obj, dofObject, index, time, vN, vN1)
            if index == 2
                obj.stateJournal(index-1).position = dofObject.qN;
                obj.stateJournal(index-1).velocity = vN;
            end
            obj.stateJournal(index).time = time;
            obj.stateJournal(index).position = dofObject.qN1;
            obj.stateJournal(index).velocity = vN1;
        end
        function out = getTime(obj, setupObject)
            out = vertcat(obj.energyJournal(:).time);
        end
        function out = getKineticEnergy(obj, setupObject)
            out = vertcat(obj.energyJournal(:).EKin);
        end
        function out = getKineticEnergyDifference(obj, setupObject)
            out = vertcat(obj.energyJournal(:).EKinDifference);
        end
        function out = getEnergy(obj, dofObject, setupObject, energyType)
            flag = false;
            out = zeros(setupObject.timeStep+1, 1);
            for ii = 1:dofObject.numberOfContinuumObjects
                continuumObject = dofObject.listContinuumObjects{ii};
                if isprop(continuumObject, 'ePot')
                    if isfield(continuumObject.ePot, energyType)
                        out = out + vertcat(continuumObject.ePot(:).(energyType));
                        flag = true;
                    end
                end
            end
            if ~flag
                out = [];
            end
        end
        function [allMaps, normMaps] = getMomentum(obj, dofObject, setupObject, momentumType, dimension)
            flag = false;
            allMaps = zeros(setupObject.timeStep+1, dimension);
            for ii = 1:dofObject.numberOfContinuumObjects
                continuumObject = dofObject.listContinuumObjects{ii};
                if isprop(continuumObject, 'momenta')
                    allMaps = allMaps + vertcat(continuumObject.momenta.(momentumType));
                    flag = true;
                    normMaps = sqrt(sum(allMaps.^2, 2));
                end
            end
            if ~flag
                allMaps = [];
                normMaps = [];
            end
        end
    end
end