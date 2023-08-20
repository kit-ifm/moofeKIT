%% startUpMoofeKIT
% run('../../startUpMoofeKIT.m');                           % Laden der Ordnerstruktur. Aufruf ist nicht notwendig, wenn das Skript schon einmal ausgeführt wurde.

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.totalTimeSteps = 5; % Anzahl der Zeitschritte anpassen (Zeitschleife)
setupObject.totalTime = 1; % Gesamtzeit der Simulation festlegen
setupObject.plotObject.flag = true; % Festlegen, dass die Ergebnisse der Simulation geplottet werden

dofObject = dofClass; % required object for dof and object handling

%% continuum object
solidObject = solidClass(dofObject);
solidObject.dimension = 2; % Dimension des zu betrachtenden Körpers angeben
solidObject.shapeFunctionObject.order = 2; % Ordnung der Formfunktionen angeben
solidObject.shapeFunctionObject.numberOfGausspoints = 4; % Anzahl der Gaußpunkte anpassen
solidObject.elementDisplacementType = 'displacement'; % Art der FE-Formulierung. Hier 'displacement'.
[solidObject.meshObject.nodes, solidObject.meshObject.edof, edofNeumann] = meshCooksMembrane(4, 4, solidObject.shapeFunctionObject.order); % Knoten und edof der Cooks-Membran beschaffen
solidObject.materialObject.name = 'HookeEVZ'; % Art des Materialgesetzes. Hier das Hookesche Gesetz im Falle eines ebenen Verzerrungszustandes.
solidObject.materialObject.rho = 0; % Dichte angeben. Falls rho=0 gewählt wird, wird rein statisch (ohne Berücksichtigung von dynamischen Effekten) gerechnet.
solidObject.materialObject.E = 100; % Elastizitaetsmodul anpassen
solidObject.materialObject.nu = 0.2; % Querkontraktionszahl festlegen

% dirchlet boundary conditions
dirichletObject = dirichletClass(dofObject);
dirichletObject.masterObject = solidObject;
dirichletObject.dimension = 2; % Dimension des zu betrachtenden Körpers angeben
dirichletObject.nodeList = find(solidObject.meshObject.nodes(:, 1) == 0);
dirichletObject.nodalDof = [1, 2];

% neumann boundary conditions
neumannObject = neumannClass(dofObject);
neumannObject.loadGeometry = 'line'; % Form des Neumann-Rand => line oder area
neumannObject.masterObject = solidObject;
neumannObject.forceVector = [0; 1]; % Komponenten des Kraftvektors anpassen
neumannObject.timeFunction = @(t) t; % Zeitfunktion, mit der die Kraft aufgetragen wird, angeben
neumannObject.meshObject.edof = edofNeumann;

%% solver
runNewton(setupObject, dofObject);

%% postprocessing
yDisplacementUpperRightNode = solidObject.qN1(end, 2) - solidObject.qR(end, 2);
fprintf('\nDisplacement: %4.3f\n', yDisplacementUpperRightNode)