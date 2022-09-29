% LSHAPEELECTROTHERMOMECHANICS Script convergence test with analytical body
% force for a static coupled electro-thermo-mechanical simulation.
%
% FORMULATION
% Mixed finite element formulation for nonlinear electro-thermo-mechanical
% processes is pursued. This is accomplished under the assumption of a
% hyperelastic, isotropic material and based on a polyconvex inspired
% internal energy function. The discretization in time is pursued with an
% energy and momentum conserving scheme.
%
%
% REFERENCE
% https://doi.org/10.1016/j.cma.2017.09.020 -> chapter 4.1
% https://doi.org/10.1016/j.cma.2021.114298
% https://doi.org/10.1007/BF00913408
%
% SEE ALSO
% LShape,
% cooksMembrane
%
% CREATOR(S)
% Marlon Franke

clear all; close all;

% Simulation Parameters
errorStudy = true;
convergenceStudy = true;
post = true;

choiceOrder = 2;

% choiceFunction = 2; %1)quadratic only for X1; 2)trigonom. in all coord. 3)cubic in all coord.
choiceFunction = 2;
% choiceElement = 'P1P0';
choiceElement = 'P2P1';
% choiceElement = 'H1H0';
% choiceElement = 'H2H1';
loadEnd = 1;
loadStepSize = 1;

% plotVar = 'stress';
% plotVar = 'phi';
plotVar = 'temp';

% Amount of Elements for Convergence-Study
All_Elem = 2:2:4;
analyticalVariablesCell = cell(1,size(All_Elem,2));
numericalVariablesCell = cell(1,size(All_Elem,2));
journalSolidElectroThermoObject = cell(1,size(All_Elem,2));

a = 1;
b = 1;
c1 = 1;
d1 = 2*(a+2*b);
dimension = 3;

% Analytical Loads (body forcy and charge density)
if choiceFunction == 3  %after Ortigosa et al. 2016
    g1 = 0.01;
    g2 = 0.02;
    g3 = 0.03;
else %e.g. Poya et al.
    g1 = 0.1/2;
    g2 = 0.2/2;
    g3 = 0.3/2;
end
[xFunction,CFunction,GFunction,cFunction,lambdaCFunction,lambdaGFunction,lambdacFunction,BFunction] = staticConvergenceAnalyticalComputations(a,b,c1,d1,g1,g2,g3,choiceFunction,dimension);
%% Loop over multiple Discretisations
for j = 1:size(All_Elem,2)
    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'solid';
    setupObject.saveObject.saveData = true;
    setupObject.totalTimeSteps = loadEnd/loadStepSize;
    setupObject.totalTime = loadEnd;
    setupObject.newtonTolerance = 1e-6;
    setupObject.plotObject.flag = false;
    setupObject.plotObject.postPlotType = 'stress';
    setupObject.plotObject.stress.component = -1;
    % setupObject.plotObject.colorBarLimits = [3200,3300];
    setupObject.plotObject.makeMovie = false;
    setupObject.integrator = 'Endpoint';

    dofObject = dofClass; % required object for dof and object handling

    %% continuum Objects
    solidObject = solidClass(dofObject);
    if  strcmpi(choiceElement,'P1P0')
        [solidObject.meshObject.nodes,solidObject.meshObject.edof,edofneumannTet] = trilinearTetrahedralBrick(1, 1, 1,All_Elem(j), All_Elem(j), All_Elem(j));
        solidObject.meshObject.nodes = solidObject.meshObject.nodes + 0.5;
        solidObject.ORDER = 1;
        solidObject.NGP = 5; %bugs for choiceFunction == 2 || 3
        solidObject.ORDER_MIXED_VAR_EL= [0 0];
    elseif strcmpi(choiceElement,'P2P1')
        [solidObject.meshObject.nodes,solidObject.meshObject.edof,edofneumannTet] = triquadraticTetrahedralBrick(All_Elem(j), All_Elem(j), All_Elem(j), 1, 1, 1, 1);
        solidObject.meshObject.nodes(:,3) = solidObject.meshObject.nodes(:,3) + 1;
        solidObject.shapeFunctionObject.order = 2;
        solidObject.shapeFunctionObject.numberOfGausspoints = 11;
    elseif strcmpi(choiceElement,'H1H0')
%         [solidElectroThermoObject.meshObject.nodes,solidElectroThermoObject.meshObject.edof] = brick3D(1, 1, 1, All_Elem(j), All_Elem(j), All_Elem(j));
        serendipity = false;
        [solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshGeneratorCube(1, 1, 1, All_Elem(j), All_Elem(j), All_Elem(j), choiceOrder, serendipity);
        solidObject.meshObject.nodes = solidObject.meshObject.nodes + 0.5;
        solidObject.shapeFunctionObject.order = 1;
        solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    elseif strcmpi(choiceElement,'H2H1')
        serendipity = true;
        [solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshGeneratorCube(1, 1, 1, All_Elem(j), All_Elem(j), All_Elem(j), choiceOrder, serendipity);
        solidObject.meshObject.nodes = solidObject.meshObject.nodes + 0.5;
        solidObject.shapeFunctionObject.order = 2;
        solidObject.shapeFunctionObject.numberOfGausspoints = 27;
    else
        error('Element-type not implemented')
    end
    solidObject.elementDisplacementType = 'mixedSC';
    solidObject.materialObject.name = 'MooneyRivlin';

    %% Material model
    solidObject.materialObject.rho = 0; % mass-density
    solidObject.materialObject.a = a;
    solidObject.materialObject.b = b;
    solidObject.materialObject.c = c1;
    solidObject.materialObject.d = d1;
    solidObject.dimension = 3;
    solidObject.mixedFEObject.condensation = true;
    solidObject.mixedFEObject.orderShapeFunction = 1;
    solidObject.numericalTangentObject.computeNumericalTangent = false;
    solidObject.numericalTangentObject.showDifferences = false;
    nodes = solidObject.meshObject.nodes;

    %% Body force
    bodyLoad = bodyForceClass(dofObject);
    bodyLoad.masterObject = solidObject;
    bodyLoad.meshObject.nodes = nodes;
    bodyLoad.meshObject.edof = solidObject.meshObject.edof;
    bodyLoad.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
    bodyLoad.shapeFunctionObject.numberOfGausspoints = solidObject.shapeFunctionObject.numberOfGausspoints;
    bodyLoad.typeOfLoad = 'deadLoad';
    bodyLoad.timeFunction = @(t) t;
    bodyLoad.loadFunction = BFunction;

    %% Dirichlet BC
    dirichletSurfaceX = dirichletClass(dofObject);
    dirichletSurfaceX.nodalDof = 1;
    dirichletSurfaceX.masterObject = solidObject;
    auessereKnoten = unique([find(nodes(:,1) == 1),find(nodes(:,1) == 0),find(nodes(:,2) == 1),find(nodes(:,2) == 0),find(nodes(:,3) == 1),find(nodes(:,3) == 0)]);
    dirichletSurfaceX.nodeList = auessereKnoten;
    if choiceFunction == 1
        dirichletSurfaceX.timeFunction = str2func('@(t,x) x+t*0.1/2*x.^2');
    elseif choiceFunction == 2
        dirichletSurfaceX.timeFunction = str2func('@(t,x) x+t*(0.1/2*sin(x))');
    elseif choiceFunction == 3
        dirichletSurfaceX.timeFunction = str2func('@(t,x) x+t*0.01*x.^3');
    elseif choiceFunction == 4
        dirichletSurfaceX.timeFunction = str2func('@(t,x) x+t*0.1/2*x.^2');
    elseif choiceFunction == 5
        dirichletSurfaceX.timeFunction = str2func('@(t,x) x');
    end

    dirichletSurfaceY = dirichletClass(dofObject);
    dirichletSurfaceY.nodalDof = 2;
    dirichletSurfaceY.masterObject = solidObject;
    dirichletSurfaceY.nodeList = auessereKnoten;
    if choiceFunction == 1
        dirichletSurfaceY.timeFunction = str2func('@(t,y) y');
    elseif choiceFunction == 2
        dirichletSurfaceY.timeFunction = str2func('@(t,y) y+t*(0.2/2*cos(y))');
    elseif choiceFunction == 3
        dirichletSurfaceY.timeFunction = str2func('@(t,y) y+t*0.02*y.^3');
    elseif choiceFunction == 4
        dirichletSurfaceY.timeFunction = str2func('@(t,y) y+t*0.2/2*y.^2');
    elseif choiceFunction == 5
        dirichletSurfaceY.timeFunction = str2func('@(t,y) y');
    end

    dirichletSurfaceZ = dirichletClass(dofObject);
    dirichletSurfaceZ.nodalDof = 3;
    dirichletSurfaceZ.masterObject = solidObject;
    dirichletSurfaceZ.nodeList = auessereKnoten;
    if choiceFunction == 1
        dirichletSurfaceZ.timeFunction = str2func('@(t,z) z');
    elseif choiceFunction == 2
        dirichletSurfaceZ.timeFunction = str2func('@(t,z) z+t*(0.3/2*(sin(z)+cos(z)))');
    elseif choiceFunction == 3
        dirichletSurfaceZ.timeFunction = str2func('@(t,z) z+t*0.03*z.^3');
    elseif choiceFunction == 4
        dirichletSurfaceZ.timeFunction = str2func('@(t,z) z+t*0.3/2*z.^2');
    elseif choiceFunction == 5
        dirichletSurfaceZ.timeFunction = str2func('@(t,z) z');
    end
    warning off
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix')
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix')
    dofObject = runNewton(setupObject,dofObject);

    if errorStudy
        %% Calculate analytic quantities for Postprocessing
        % motion
        analyticalVariables.x = zeros(size(nodes(:,1:3)));
        % strain measures
        analyticalVariables.C = zeros(size(nodes,1),6);
        analyticalVariables.G = zeros(size(nodes,1),6);
        analyticalVariables.c = zeros(size(nodes,1),1);
        % stress type LM
        analyticalVariables.lambdaC = zeros(size(nodes,1),6);
        analyticalVariables.lambdaG = zeros(size(nodes,1),6);
        analyticalVariables.lambdac = zeros(size(nodes,1),1);
        for i=1:size(analyticalVariables.x,1)
            analyticalVariables.x(i,:)          = xFunction(nodes(i,1), nodes(i,2),nodes(i,3));
            CTensor       = CFunction(nodes(i,1), nodes(i,2), nodes(i,3));
            analyticalVariables.C(i,:) = [CTensor(1,1) CTensor(2,2) CTensor(3,3) CTensor(1,2) CTensor(2,3) CTensor(1,3)];
            GTensor       = GFunction(nodes(i,1), nodes(i,2), nodes(i,3));
            analyticalVariables.G(i,:) = [GTensor(1,1) GTensor(2,2) GTensor(3,3) GTensor(1,2) GTensor(2,3) GTensor(1,3)];
            analyticalVariables.c(i,:) = cFunction(nodes(i,1), nodes(i,2), nodes(i,3));
            lambdaCTensor  = lambdaCFunction(nodes(i,1), nodes(i,2), nodes(i,3));
            analyticalVariables.lambdaC(i,:)  = [lambdaCTensor(1,1) lambdaCTensor(2,2) lambdaCTensor(3,3) lambdaCTensor(1,2) lambdaCTensor(2,3) lambdaCTensor(1,3)];
            lambdaGTensor  = lambdaGFunction(nodes(i,1), nodes(i,2), nodes(i,3));
            analyticalVariables.lambdaG(i,:)  = [lambdaGTensor(1,1) lambdaGTensor(2,2) lambdaGTensor(3,3) lambdaGTensor(1,2) lambdaGTensor(2,3) lambdaGTensor(1,3)];
            analyticalVariables.lambdac(i,:)  = lambdacFunction(nodes(i,1), nodes(i,2), nodes(i,3));
        end
        analyticalVariablesCell{j}.q      = analyticalVariables.x;
        analyticalVariablesCell{j}.C      = analyticalVariables.C;
        analyticalVariablesCell{j}.G      = analyticalVariables.G;
        analyticalVariablesCell{j}.c      = analyticalVariables.c;
        analyticalVariablesCell{j}.lambdaC = analyticalVariables.lambdaC;
        analyticalVariablesCell{j}.lambdaG = analyticalVariables.lambdaG;
        analyticalVariablesCell{j}.lambdac = analyticalVariables.lambdac;
        journalSolidElectroThermoObject{j} = solidObject;
        Elem{1}.disp(j,:) = solidObject.qN1(find(nodes(:,1) == 0.5 & nodes(:,2) == 0.5 & nodes(:,3) == 0.5),1:3) - 0.5;
    end
end

if post == true
    staticConvergencePost;
end

if convergenceStudy
    plotcolors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980  ];  % b/r
    plot(All_Elem,Elem{1}.disp(:,2),'Color',plotcolors(1,:),'linewidth',1.2);
    hold on
    plot(All_Elem,Elem{2}.disp(:,2),'Color',plotcolors(2,:),'linewidth',1.2)
    hold off
    xlabel('n_{el} per side')
    ylabel('u_2(0.5,0.5,0.5)')
    legend('mmhwSC','mmhwCascadeSC')
    grid on
end
