% run solid examples


% add folders to workspace
% run('../startUpMoofeKIT.m')

% start parpool
% parpool(24)

disp('=======================================');
disp('Now testing: GFE scripts');
disp('=======================================');

% turn off figures
set(0,'DefaultFigureVisible','off');


%% Test 1: check displacement of 'Aufgabe13a' from lecture GFE
run('../scripts/pre/GFE/Aufgabe13a.m')
assert(displacementLastNodeArray(1, 1) - (-11.9047619) <= 1e-7, 'Error in GFE/Aufgabe13a');
assert(displacementLastNodeArray(2, 1) - (-0.1543e-3) <= 1e-7, 'Error in GFE/Aufgabe13a');
% assert(displacementLastNodeArray(1, 2) - (-0.301759) <= 1e-6, 'Error in GFE/Aufgabe13a');
% assert(displacementLastNodeArray(2, 2) - (-1.6520e-04) <= 1e-8, 'Error in GFE/Aufgabe13a');

% turn on figures
set(0,'DefaultFigureVisible','on');
