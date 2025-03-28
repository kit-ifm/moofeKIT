%% startUpMoofeKIT
% run('../../startUpMoofeKIT.m');                           % Laden der Ordnerstruktur. Aufruf ist nicht notwendig, wenn das Skript schon einmal ausgeführt wurde.

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.totalTimeSteps = 1; % Anzahl der Zeitschritte anpassen (Zeitschleife)
setupObject.totalTime = 1; % Gesamtzeit der Simulation festlegen
setupObject.plotObject.flag = false; % Festlegen, dass die Ergebnisse der Simulation geplottet werden
setupObject.plotObject.view = 2; % Plotten als 2D-Objekt

dofObject = dofClass; % required object for dof and object handling

%% continuum object
solidObject = solidClass(dofObject);
solidObject.dimension = 2; % Dimension des zu betrachtenden Körpers angeben
solidObject.shapeFunctionObject.order = 1; % Ordnung der Formfunktionen angeben
solidObject.shapeFunctionObject.numberOfGausspoints = 4; % Anzahl der Gaußpunkte anpassen
solidObject.elementDisplacementType = 'displacement'; % Art der FE-Formulierung. Hier 'displacement'.
solidObject.elementNameAdditionalSpecification = 'PetrovGalerkinXieEtAl2016TQ4'; % Konkretisierung Elementname
[solidObject.meshObject.nodes, solidObject.meshObject.edof, edofNeumann] = meshRectangle(2, 2, 1, 1, solidObject.shapeFunctionObject.order, false); % Knoten und edof der Cooks-Membran beschaffen
solidObject.meshObject.nodes(:, 1) = solidObject.meshObject.nodes(:, 1) + 1;
solidObject.meshObject.nodes(:, 2) = solidObject.meshObject.nodes(:, 2) + 1;
solidObject.materialObject.name = 'HookeESZ'; % Art des Materialgesetzes. Hier das Hookesche Gesetz im Falle eines ebenen Verzerrungszustandes.
solidObject.materialObject.rho = 0; % Dichte angeben. Falls rho=0 gewählt wird, wird rein statisch (ohne Berücksichtigung von dynamischen Effekten) gerechnet.
solidObject.materialObject.E = 1; % Elastizitaetsmodul anpassen
solidObject.materialObject.nu = 1/3; % Querkontraktionszahl festlegen

% dirchlet boundary conditions
dirichletObject = dirichletClass(dofObject);
dirichletObject.masterObject = solidObject;
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
dirichletObject.nodalDof = [1, 2];

% neumann boundary conditions
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'line'; % Form des Neumann-Rand => line oder area
neumannObject.masterObject = solidObject;
neumannObject.loadVector = [0; 1]/16; % Komponenten des Lastvektors anpassen
neumannObject.meshObject.edof = edofNeumann.SX2;

%% solver
runNewton(setupObject, dofObject);

%% postprocessing
nodeToMeasure = find(abs(solidObject.meshObject.nodes(:, 1)-2) <= 1e-8 & abs(solidObject.meshObject.nodes(:, 2)-2) <= 1e-8);
yDisplacementPointC = solidObject.qN1(nodeToMeasure, 2) - solidObject.qN(nodeToMeasure, 2);
fprintf('\nNormalized Displacement: %4.3f\n', yDisplacementPointC/23.96);