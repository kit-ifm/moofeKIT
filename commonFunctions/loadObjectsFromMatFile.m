function [setupObjectCell, dofObjectCell] = loadObjectsFromMatFile(filename, timeSteps)
%LOADOBJECTSFROMMATFILE Loads setupObjects and dofObjects from .mat file
%   This function loads the moofeKIT objects 'setupObject' and
%   'dofObject' for all given time steps and returns them in cell arrays.
%
%   CALL
%   [setupObjectCell, dofObjectCell] = loadObjectsFromMatFile(filename, timeSteps)
%   setupObjectCell:    cell containing all requested setupObjects
%   dofObjectCell:      cell containing all requested dofObjects
%   filename:           string (e.g. 'fileXYZ')
%   timeSteps:          vector with required time steps or 'all' or 'first' or 'last'
%
%   CREATOR(S)
%   Felix Zaehringer

% check input
assert(ischar(filename), 'filename must be of type string!');
assert((isvector(timeSteps) && isnumeric(timeSteps)) || any(strcmpi(timeSteps, {'all', 'first', 'last'})), 'timeSteps must be a vector or a string of the set {"all", "first", "last"}!');
generatedFilesFolderInfo = what('generatedFiles');
filename = strcat(generatedFilesFolderInfo.path, filesep, 'matFiles', filesep, filename);
if ~strcmp(filename(end-3:end), '.mat')
    filename = strcat(filename, '.mat');
end
assert(isfile(filename), 'The requested .mat file could not be found!');

% find out which time steps are stored in the .mat file
timeStepsInFile = sort(str2double(unique(extract(who('-file', filename), digitsPattern))));

% determine which variables / objects should be loaded
if isvector(timeSteps) && isnumeric(timeSteps)
    requestedTimeSteps = timeSteps;
elseif strcmpi(timeSteps, 'all')
    requestedTimeSteps = timeStepsInFile;
elseif strcmpi(timeSteps, 'first')
    requestedTimeSteps = timeStepsInFile(1);
elseif strcmpi(timeSteps, 'last')
    requestedTimeSteps = timeStepsInFile(end);
end
numberOfTimeStepsRequested = length(requestedTimeSteps);

variablesToLoad = cell(2*numberOfTimeStepsRequested, 1);
for i = 1:numberOfTimeStepsRequested
    if ismember(requestedTimeSteps(i), timeStepsInFile)
        index = 2 * i - 1;
        variablesToLoad{index} = strcat('setupObject', num2str(requestedTimeSteps(i)));
        variablesToLoad{index+1} = strcat('dofObject', num2str(requestedTimeSteps(i)));
    else
        error(filename+" does not contain setupObject or dofObject for time step "+num2str(requestedTimeSteps(i))+"!");
    end
end

% load variables / objects
load(filename, variablesToLoad{:});

% save loaded variables / objects in cell array
setupObjectCell = cell(numberOfTimeStepsRequested, 1);
dofObjectCell = cell(numberOfTimeStepsRequested, 1);
for i = 1:numberOfTimeStepsRequested
    eval(strcat('setupObjectCell{i} = setupObject', num2str(requestedTimeSteps(i)), ';'));
    eval(strcat('dofObjectCell{i} = dofObject', num2str(requestedTimeSteps(i)), ';'));
end

end
