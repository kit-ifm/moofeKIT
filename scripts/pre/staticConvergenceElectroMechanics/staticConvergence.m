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

choiceFunction = 2; %1)quadratic only for X1; 2)trigonom. in all coord. 3)cubic in all coord.
% choiceFunction = 3;
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
journalSolidElectroObject = cell(1,size(All_Elem,2));

a = 1;
b = 1;
c = 1;
d = 2*(a+2*b);
eps_0 = 8.854*10^(-12);
eps_r = 4;
dimension = 3;
phi0 = 10^2;

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
[xFunction,phiFunction,DFunction,BSymbolic,rho0Symbolic,CFunction,GFunction,cFunction,lambdaCFunction,lambdaGFunction,lambdacFunction] = staticConvergenceAnalyticalComputations(a,b,c,d,eps_0,eps_r,g1,g2,g3,phi0,choiceFunction,dimension);
for j = 1:size(All_Elem,2)
    %% setup (mandatory: setup and dofs)
    setupObject = setupClass;
    setupObject.saveObject.fileName = 'solidElectroThermo';
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
    solidElectroObject = solidElectroClass(dofObject);
    if  strcmpi(choiceElement,'P1P0')
        [solidElectroObject.meshObject.nodes,solidElectroObject.meshObject.edof,edofneumannTet] = trilinearTetrahedralBrick(1, 1, 1,All_Elem(j), All_Elem(j), All_Elem(j));
        solidElectroObject.meshObject.nodes = solidElectroObject.meshObject.nodes + 0.5;
        solidElectroObject.ORDER = 1;
        solidElectroObject.NGP = 5; %bugs for choiceFunction == 2 || 3
        solidElectroObject.ORDER_MIXED_VAR_EL= [0 0];
    elseif strcmpi(choiceElement,'P2P1')
        [solidElectroObject.meshObject.nodes,solidElectroObject.meshObject.edof,edofneumannTet] = triquadraticTetrahedralBrick(All_Elem(j), All_Elem(j), All_Elem(j), 1, 1, 1, 1);
        solidElectroObject.meshObject.nodes(:,3) = solidElectroObject.meshObject.nodes(:,3) + 1;
        solidElectroObject.shapeFunctionObject.order = 2;
        solidElectroObject.shapeFunctionObject.numberOfGausspoints = 11;
        solidElectroObject.mixedFEObject.numberOfNodes = 4;
    elseif strcmpi(choiceElement,'H1H0')
%         [solidElectroThermoObject.meshObject.nodes,solidElectroThermoObject.meshObject.edof] = brick3D(1, 1, 1, All_Elem(j), All_Elem(j), All_Elem(j));
        serendipity = false;
        [solidElectroObject.meshObject.nodes,solidElectroObject.meshObject.edof] = meshGeneratorCube(1, 1, 1, All_Elem(j), All_Elem(j), All_Elem(j), choiceOrder, serendipity);
        solidElectroObject.meshObject.nodes = solidElectroObject.meshObject.nodes + 0.5;
        solidElectroObject.shapeFunctionObject.order = 1;
        solidElectroObject.shapeFunctionObject.numberOfGausspoints = 8;
    elseif strcmpi(choiceElement,'H2H1')
%         serendipity = true;
%         [solidElectroObject.meshObject.nodes,solidElectroObject.meshObject.edof] = meshGeneratorCube(1, 1, 1, All_Elem(j), All_Elem(j), All_Elem(j), choiceOrder, serendipity);
%         solidElectroObject.meshObject.nodes = solidElectroObject.meshObject.nodes + 0.5;
        [solidElectroObject.meshObject.nodes,solidElectroObject.meshObject.edof] = triquadraticBrick(All_Elem(j), All_Elem(j), All_Elem(j), 1, 1, 1, 1);
        solidElectroObject.shapeFunctionObject.order = 2;
        solidElectroObject.shapeFunctionObject.numberOfGausspoints = 27;
    else
        error('Element-type not implemented')
    end
    initialElectricalPotentialField = zeros(size(solidElectroObject.meshObject.nodes,1),1);
    solidElectroObject.meshObject.nodes = [solidElectroObject.meshObject.nodes, initialElectricalPotentialField];
%     solidElectroObject.elementDisplacementType = 'mixedD_SC';
%     solidElectroObject.materialObject.name = 'MooneyRivlin';
    solidElectroObject.elementDisplacementType = 'mixedSC';
    solidElectroObject.materialObject.name = 'MooneyRivlin';

    %% Material model
    solidElectroObject.materialObject.rho = 0; % mass-density
    solidElectroObject.materialObject.a = a;
    solidElectroObject.materialObject.b = b;
    solidElectroObject.materialObject.c = c;
    solidElectroObject.materialObject.d = d;
    solidElectroObject.materialObject.e0 = eps_0; % epsilon_0
    solidElectroObject.materialObject.e1 = eps_r; % epsilon_r
    solidElectroObject.materialObject.e2 = 0;
    solidElectroObject.materialObject.ee = 0;
    solidElectroObject.materialObject.rhoSource = rho0Symbolic;
    solidElectroObject.materialObject.timeFunctionRhoSource = str2func('@(t) t');
    solidElectroObject.dimension = 3;
    solidElectroObject.mixedFEObject.condensation = true;
    solidElectroObject.mixedFEObject.orderShapeFunction = 1;
    solidElectroObject.numericalTangentObject.computeNumericalTangent = false;
    solidElectroObject.numericalTangentObject.showDifferences = false;
    nodes = solidElectroObject.meshObject.nodes;

    %% Body force
    bodyLoad = bodyForceClass(dofObject);
    bodyLoad.masterObject = solidElectroObject;
    bodyLoad.meshObject.nodes = nodes;
    bodyLoad.meshObject.edof = solidElectroObject.meshObject.edof;
    bodyLoad.shapeFunctionObject.order = solidElectroObject.shapeFunctionObject.order;
    bodyLoad.shapeFunctionObject.numberOfGausspoints = solidElectroObject.shapeFunctionObject.numberOfGausspoints;
    bodyLoad.typeOfLoad = 'deadLoad';
    bodyLoad.timeFunction = @(t) t;
    bodyLoad.loadFunction = matlabFunction(formula(BSymbolic));

    %% Dirichlet BC
    dirichletSurfaceX = dirichletClass(dofObject);
    dirichletSurfaceX.nodalDof = 1;
    dirichletSurfaceX.masterObject = solidElectroObject;
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
    dirichletSurfaceY.masterObject = solidElectroObject;
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
    dirichletSurfaceZ.masterObject = solidElectroObject;
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

    dirichletElectricalPotential = dirichletClass(dofObject);
    dirichletElectricalPotential.functionDof = 1;%'addfield1X';
    dirichletElectricalPotential.nodalDof = solidElectroObject.dimension + 1;
    dirichletElectricalPotential.masterObject = solidElectroObject;
    dirichletElectricalPotential.nodeList = auessereKnoten;
    if choiceFunction == 1
        dirichletElectricalPotential.timeFunction = str2func('@(t,X1) t*10^4 * sin(X1)');
    elseif choiceFunction == 2
        dirichletElectricalPotential.timeFunction = str2func('@(t,X1) t*10^2 * sin(X1)');
    elseif choiceFunction == 3
        dirichletElectricalPotential.timeFunction = str2func('@(t,X1) t*10^2 * X1.^3');
    elseif choiceFunction == 4
        dirichletElectricalPotential.timeFunction = str2func('@(t,X1) t*10^2 * X1');
    elseif choiceFunction == 5
        dirichletElectricalPotential.timeFunction = str2func('@(t,X1) t*0 * X1');
    end

    %% Test
    % xanaly
    if choiceFunction == 2
        for indexX = 0:0.5:1
            if abs(double(phiFunction(indexX,0,0) - dirichletElectricalPotential.timeFunction(1,indexX))) >= 1e-12
                error('Analytical formulations are divergent')
            end
            for indexY = 0:0.5:1
                for indexZ = 0:0.5:1
                    if abs(double(xFunction(indexX,indexY,indexZ)-[dirichletSurfaceX.timeFunction(1,indexX);dirichletSurfaceY.timeFunction(1,indexY);dirichletSurfaceZ.timeFunction(1,indexZ)])) >= 1e-12
                        error('Analytical formulations are divergent')
                    end
                end
            end
        end
    end
    if choiceFunction == 4
        for indexX = 0:0.5:1
            if abs(double(phiFunction(indexX,0,0) - dirichletElectricalPotential.timeFunction(1,indexX))) >= 1e-12
                error('Analytical formulations are divergent')
            end
            for indexY = 0:0.5:1
                for indexZ = 0:0.5:1
                    if abs(double(xFunction(indexX,indexY,indexZ)-[dirichletSurfaceX.timeFunction(1,indexX);dirichletSurfaceY.timeFunction(1,indexY);dirichletSurfaceZ.timeFunction(1,indexZ)])) >= 1e-12
                        error('Analytical formulations are divergent')
                    end
                end
            end
        end
    end

    warning off
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix')
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix')
    dofObject = runNewton(setupObject,dofObject);

    if errorStudy
        %% Calculate analytic quantities for Postprocessing
        % motion
        analyticalVariables.x = zeros(size(nodes(:,1:3)));
        % electric potential
        analyticalVariables.phi = zeros(size(nodes(:,1)));
        % electrical Displacement
        analyticalVariables.D = zeros(size(nodes(:,1:3)));
        if strcmpi(solidElectroObject.elementDisplacementType,'mixedSC')
            % strain measures
            analyticalVariables.C = zeros(size(nodes,1),6);
            analyticalVariables.G = zeros(size(nodes,1),6);
            analyticalVariables.c = zeros(size(nodes,1),1);
            % stress type LM
            analyticalVariables.lambdaC = zeros(size(nodes,1),6);
            analyticalVariables.lambdaG = zeros(size(nodes,1),6);
            analyticalVariables.lambdac = zeros(size(nodes,1),1);
        end
        for i=1:size(analyticalVariables.x,1)
            analyticalVariables.x(i,:)          = xFunction(nodes(i,1), nodes(i,2),nodes(i,3));
            analyticalVariables.phi(i,:)        = phiFunction(nodes(i,1),nodes(i,2),nodes(i,3));
            analyticalVariables.D(i,:)          = DFunction(nodes(i,1),nodes(i,2),nodes(i,3));
            if strcmpi(solidElectroObject.elementDisplacementType,'mixedSC')
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
        end
        analyticalVariablesCell{j}.q      = [analyticalVariables.x,analyticalVariables.phi];
        analyticalVariablesCell{j}.D      = analyticalVariables.D;
        if strcmpi(solidElectroObject.elementDisplacementType,'mixedSC')
            analyticalVariablesCell{j}.C      = analyticalVariables.C;
            analyticalVariablesCell{j}.G      = analyticalVariables.G;
            analyticalVariablesCell{j}.c      = analyticalVariables.c;
            analyticalVariablesCell{j}.lambdaC = analyticalVariables.lambdaC;
            analyticalVariablesCell{j}.lambdaG = analyticalVariables.lambdaG;
            analyticalVariablesCell{j}.lambdac = analyticalVariables.lambdac;
        end
        journalSolidElectroObject{j} = solidElectroObject;
        Elem{1}.disp(j,:) = solidElectroObject.qN1(find(nodes(:,1) == 0.5 & nodes(:,2) == 0.5 & nodes(:,3) == 0.5),1:3) - 0.5;
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
