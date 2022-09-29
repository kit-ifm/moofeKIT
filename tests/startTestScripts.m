% startScript for tests

% add folders to workspace
run('../startUpMoofeKIT.m')

% start parpool
% parpool(24)

result = runtests('testMeSolidElements.m');