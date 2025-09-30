% run solid examples


% add folders to workspace
% run('../startUpMoofeKIT.m')

% start parpool
% parpool(24)

% turn off figures
set(0,'DefaultFigureVisible','off')

%% Test 1: check the total energy of a mixed thermo-electro-elastodynamic example
run('../scripts/pre/unitTest.m')
assert(lastTotalEnergyValue - (-1.334785464144597e+04) < 1e-15, 'Total energy in unitTest.m does not equal the specified total energy');

%% Test 2: check 'Lager' load of 'Aufgabe02' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe02.m')
assert((sum(lager)>9.99 & sum(lager) < 10.01),'sum of lager stress should give a load of 10')

%% Test 3: check y-displacement of upper right node of 'Aufgabe03' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe03.m')
assert(yDisplacementUpperRightNode > 15.010 & yDisplacementUpperRightNode < 15.012,'y-displacement of upper right node does not correspond to the expected value')

%% Test 4: check y-displacement of upper right node of 'Aufgabe04' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe04.m')
assert(yDisplacementUpperRightNode > 0.658 & yDisplacementUpperRightNode < 0.66,'y-displacement of upper right node does not correspond to the expected value')

%% Test 5: check y-displacement of upper right node of 'Aufgabe05' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe05.m')
assert(yDisplacementUpperRightNode > 17.2739 & yDisplacementUpperRightNode < 17.2741,'y-displacement of upper right node does not correspond to the expected value')

%% Test 6: check y-displacement of upper right node of 'Aufgabe06' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe06.m')
assert(yDisplacementUpperRightNode > 15.7034 & yDisplacementUpperRightNode < 15.7036,'y-displacement of upper right node does not correspond to the expected value')

%% Test 7: check y-displacement of upper right node of 'Aufgabe07' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe07.m')
assert(min(yDisplacementUpperRightNode) > 5.629 & min(yDisplacementUpperRightNode) < 5.631,'y-displacement of upper right node does not correspond to the expected value')

%% Test 8: check y-displacement of upper right node of 'Aufgabe08' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe08.m')
assert(yDisplacementUpperRightNode > 14.9767 & yDisplacementUpperRightNode < 14.9769,'y-displacement of upper right node does not correspond to the expected value')

%% Test 9: check y-displacement of upper right node of 'Aufgabe09' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe09.m')
assert(yDisplacementUpperRightNode > 14.9167 & yDisplacementUpperRightNode < 14.9169,'y-displacement of upper right node does not correspond to the expected value')

%% Test 10: check y-displacement of upper right node of 'Aufgabe10' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe10.m')
assert(yDisplacementUpperRightNode > 14.6723 & yDisplacementUpperRightNode < 14.6725,'y-displacement of upper right node does not correspond to the expected value')

%% Test 11: check y-displacement of upper right node of 'Aufgabe11' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe11.m')
assert(yDisplacementUpperRightNode > 14.731 & yDisplacementUpperRightNode < 14.733,'y-displacement of upper right node does not correspond to the expected value')

%% Test 12: check y-displacement of upper right node of 'Aufgabe12' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe12.m')
assert(yDisplacementUpperRightNode > 12.061 & yDisplacementUpperRightNode < 12.063,'y-displacement of upper right node does not correspond to the expected value')

%% Test 13a: check y-displacement of upper right node of 'Aufgabe13a' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe13a.m')
assert(yDisplacementUpperRightNode > -0.236 & yDisplacementUpperRightNode < -0.234,'y-displacement of upper right node does not correspond to the expected value')

%% Test 13b: check y-displacement of upper right node of 'Aufgabe13b' from lecture FEidF
run('../scripts/pre/FEidF/Aufgabe13b.m')
assert(yDisplacementUpperRightNode > 0.0019 & yDisplacementUpperRightNode < 0.0021,'y-displacement of upper right node does not correspond to the expected value')

% turn on figures
set(0,'DefaultFigureVisible','on')
