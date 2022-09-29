classdef mapVoigtClass < matlab.mixin.Copyable
    properties (SetAccess = private)
        mapVoigt
        reMapVoigt
        mapType % 'symmetric' or 'unsymmetric'
        dimension
    end
    methods
        function selectMapVoigt(obj,dimension,mapType)
% function for Voigt mapping for 2 and 3 dimensional matrices with indices
% [ 11 12 13;   equals index: [ 1 4 7;
%   21 22 23;                   2 5 8;
%   31 32 33]                   3 6 9]
            obj.mapType = mapType;
            obj.dimension = dimension;
            if strcmpi(obj.mapType,'symmetric')
                if dimension == 1
                    obj.mapVoigt = [1];
                elseif dimension == 2
                    % chose this order: 11,22,12; corresponds to 
                    obj.mapVoigt = [1,5,4]';
                elseif dimension == 3
                    % chose this order: 11,22,33,12,23,13; corresponds to 
                    obj.mapVoigt = [1;5;9;4;8;7];
                    obj.reMapVoigt = [];
                else 
                    error('dimension not implemented');
                end
            elseif strcmpi(obj.mapType,'unsymmetric')
                %voigt map
                if dimension == 1
                    obj.mapVoigt = [1];
                elseif dimension == 2
                    % chose this order: 11,22,12,21; corresponds to 
                    obj.mapVoigt = [1,5,4,2]';
                elseif dimension == 3
                    % chose this order: 11,22,33,12,23,13,21,32,31; corresponds to 
                    obj.mapVoigt = [1;5;9;4;8;7;2;6;3];
                else 
                    error('dimension not implemented');                    
                end
            else 
                error('Map Voigt type not implemented');
            end
        end
        function out = get.mapVoigt(obj)
            out = obj.mapVoigt;
        end
    end
end