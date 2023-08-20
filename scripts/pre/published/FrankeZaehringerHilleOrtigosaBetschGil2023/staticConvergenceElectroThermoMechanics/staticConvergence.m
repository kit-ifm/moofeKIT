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

% choiceFunction = 2; %1)quadratic only for X1; 2)trigonom. in all coord. 3)cubic in all coord.
choiceFunction = 3;
% choiceElement = {'P1P0','P2P1','H1H0','H2H1'};
choiceElement = {'P2P1','H1H0','H2H1'};
% choiceElement = {'P2P1','H2H1'};
% choiceElement = {'H1H0'};
% choiceElement = {'P1P0'};
loadEnd = 1;
loadStepSize = 0.1;
% loadStepSize = 1;

% Amount of Elements for Convergence-Study
All_Elem = 4:2:8;
analyticalVariablesCell = cell(numel(choiceElement),size(All_Elem,2));
solidElectroThermoObjectCell = cell(numel(choiceElement),size(All_Elem,2));

% a = 1;
% b = 1;
% c1 = 1;
% d1 = 2*(a+2*b);
% c2 = 1;
% d2 = 0;
% eps_0 = 8.854*10^(-12);
% eps_r = 4;
% kappa = 1;
% beta = 2.233*10^(-4);                         % coupling parameter
% k0 = 0.1;                                      % thermal conductivity
% wk = 0;
a = 25000;
b = 50000;
c1 = 500000;
d1 = 2*(a+2*b);
c2 = 5209;
d2 = 0;
eps_0 = 8.854*10^(-12);
eps_r = 4;
kappa = 1500;
beta = 2.233*10^(-4);                         % coupling parameter
k0 = 0.23;                                      % thermal conductivity
wk = 0;

dimension = 3;
phi0 = 10^2;
thetaR = 293.15;                            % Reference Temperature
theta0 = 293.15;                            % Initial Temperature
thetaRb = 10;                              % Reference Temperature

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
[xFunction,phiFunction,thetaFunction,DFunction,CFunction,GFunction,cFunction,lambdaCFunction,lambdaGFunction,lambdacFunction,BSymbolic,rho0Symbolic,RSymbolic] = staticConvergenceAnalyticalComputations(a,b,c1,d1,c2,d2,eps_0,eps_r,kappa,beta,k0,g1,g2,g3,phi0,thetaR,thetaRb,choiceFunction,dimension);
%% Loop over multiple Discretisations
for ii = 1:numel(choiceElement)
    for j = 1:size(All_Elem,2)
        %% setup (mandatory: setup and dofs)
        setupObject = setupClass;
        setupObject.totalTimeSteps = loadEnd/loadStepSize;
        setupObject.saveObject.fileName = strcat(mfilename,choiceElement{ii});
        setupObject.saveObject.saveData = true;
        setupObject.saveObject.saveForTimeSteps = setupObject.totalTimeSteps;
        setupObject.saveObject.exportLogFile = true;
        setupObject.totalTime = loadEnd;
        setupObject.newton.tolerance = 1e-6;
        setupObject.plotObject.flag = false;
        setupObject.plotObject.postPlotType = 'stress';
        setupObject.plotObject.stress.component = -1;
        % setupObject.plotObject.colorBarLimits = [3200,3300];
        setupObject.plotObject.makeMovie = false;
        setupObject.integrator = 'Endpoint';

        dofObject = dofClass; % required object for dof and object handling

        %% continuum Objects
        solidElectroThermoObject = solidElectroThermoClass(dofObject);
        if  strcmpi(choiceElement{ii},'P1P0')
            choiceOrder = 1;
            choiceOrderMixed = 'sameOrder';
            choiceOrderMixedData = 0;
            [solidElectroThermoObject.meshObject.nodes,solidElectroThermoObject.meshObject.edof,edofneumannTet] = trilinearTetrahedralBrick(1, 1, 1,All_Elem(j), All_Elem(j), All_Elem(j));
            solidElectroThermoObject.meshObject.nodes = solidElectroThermoObject.meshObject.nodes + 0.5;
%             solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 5; % bugs in esra-code for choiceFunction == 2 || 3 % TODO: number of Gausspoints (5) not implemented!
%             solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 11;
            solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 5;
            solidElectroThermoObject.mixedFEObject.numberOfNodes = 1;
        elseif strcmpi(choiceElement{ii},'P2P1')
            choiceOrder = 2;
            choiceOrderMixed = 'sameOrder';
            choiceOrderMixedData = 1;
            [solidElectroThermoObject.meshObject.nodes,solidElectroThermoObject.meshObject.edof,edofneumannTet] = triquadraticTetrahedralBrick(All_Elem(j), All_Elem(j), All_Elem(j), 1, 1, 1, 1);
            solidElectroThermoObject.meshObject.nodes(:,3) = solidElectroThermoObject.meshObject.nodes(:,3) + 1;
            solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 11;
            solidElectroThermoObject.mixedFEObject.numberOfNodes = 4;
        elseif strcmpi(choiceElement{ii},'H1H0')
            choiceOrder = 1;
            serendipity = false;
            choiceOrderMixed = 'sameOrder';
            choiceOrderMixedData = 0;
            [solidElectroThermoObject.meshObject.nodes,solidElectroThermoObject.meshObject.edof] = meshGeneratorCube(1, 1, 1, All_Elem(j), All_Elem(j), All_Elem(j), choiceOrder, serendipity);
            solidElectroThermoObject.meshObject.nodes = solidElectroThermoObject.meshObject.nodes + 0.5;
            solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 8;
            solidElectroThermoObject.mixedFEObject.numberOfNodes = 1;
        elseif strcmpi(choiceElement{ii},'H2H1')
            choiceOrder = 2;
            serendipity = true;
%             choiceOrderMixed = 'detailedOrder';
%             choiceOrderMixedData = struct('D',[1 1],...
%                                       'C',[1 1],...
%                                       'G',[1 1],...
%                                       'c',[1 1],...
%                                       'LambdaC',[1 1],...
%                                       'LabmdaG',[1 1],...
%                                       'Lambdac',[1 1]);
%             choiceOrderMixedData = struct('D',2*[1 1],...
%                                       'C',2*[1 1],...
%                                       'G',2*[1 1],...
%                                       'c',2*[1 1],...
%                                       'LambdaC',2*[1 1],...
%                                       'LabmdaG',2*[1 1],...
%                                       'Lambdac',2*[1 1]);
%             choiceOrderMixedData = struct('D',[1 1]);
            choiceOrderMixed = 'sameOrder';
            choiceOrderMixedData = 1;
%             choiceOrderMixedData = 2;
            [solidElectroThermoObject.meshObject.nodes,solidElectroThermoObject.meshObject.edof] = meshGeneratorCube(1, 1, 1, All_Elem(j), All_Elem(j), All_Elem(j), choiceOrder, serendipity);
            solidElectroThermoObject.meshObject.nodes = solidElectroThermoObject.meshObject.nodes + 0.5;
            solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints = 27;
        else
            error('Element-type not implemented')
        end
        initialThermalField = theta0*ones(size(solidElectroThermoObject.meshObject.nodes,1),1);
        initialElectricalPotentialField = zeros(size(solidElectroThermoObject.meshObject.nodes,1),1);
        solidElectroThermoObject.meshObject.nodes = [solidElectroThermoObject.meshObject.nodes, initialElectricalPotentialField, initialThermalField];
%             solidElectroThermoObject.elementDisplacementType = 'mixedD_SC';
%             solidElectroThermoObject.materialObject.name = 'MooneyRivlin';
        solidElectroThermoObject.elementDisplacementType = 'mixedSC';
        solidElectroThermoObject.materialObject.name = 'MooneyRivlinFullCoupled';

        %% Material model
        solidElectroThermoObject.materialObject.rho = 0; % mass-density
        solidElectroThermoObject.materialObject.a = a;
        solidElectroThermoObject.materialObject.b = b;
        solidElectroThermoObject.materialObject.c1 = c1;
        solidElectroThermoObject.materialObject.c2 = c2;
        solidElectroThermoObject.materialObject.d1 = d1;
        solidElectroThermoObject.materialObject.d2 = d2;
        solidElectroThermoObject.materialObject.e0 = eps_0; % epsilon_0
        solidElectroThermoObject.materialObject.e1 = eps_r; % epsilon_r
        solidElectroThermoObject.materialObject.e2 = 0;
        solidElectroThermoObject.materialObject.ee = 0;
        solidElectroThermoObject.materialObject.kappa = kappa; % heat capacity
        solidElectroThermoObject.materialObject.strongcomp = false;
        solidElectroThermoObject.materialObject.beta = beta; % coupling parameter
        solidElectroThermoObject.materialObject.thetaR = thetaR;
        solidElectroThermoObject.materialObject.k0 = k0; % thermal conductivity
        solidElectroThermoObject.materialObject.wk = wk;
        solidElectroThermoObject.materialObject.rhoSource = rho0Symbolic;
        solidElectroThermoObject.materialObject.RSource = RSymbolic;
        timeFktRhoSource = str2func('@(t) t');
        solidElectroThermoObject.materialObject.timeFunctionRhoSource = timeFktRhoSource;
        timeFunctionRSource = str2func('@(t) t');
        solidElectroThermoObject.materialObject.timeFunctionRSource = timeFunctionRSource;
        solidElectroThermoObject.dimension = 3;
        solidElectroThermoObject.mixedFEObject.condensation = true;
        solidElectroThermoObject.mixedFEObject.typeShapeFunction = choiceOrderMixed;
        solidElectroThermoObject.mixedFEObject.typeShapeFunctionData = choiceOrderMixedData;
        solidElectroThermoObject.numericalTangentObject.computeNumericalTangent = false;
        solidElectroThermoObject.numericalTangentObject.showDifferences = false;
        nodes = solidElectroThermoObject.meshObject.nodes;

        %% Body force
        bodyLoad = bodyForceClass(dofObject);
        bodyLoad.masterObject = solidElectroThermoObject;
        bodyLoad.meshObject.nodes = nodes;
        bodyLoad.meshObject.edof = solidElectroThermoObject.meshObject.edof;
        bodyLoad.shapeFunctionObject.numberOfGausspoints = solidElectroThermoObject.shapeFunctionObject.numberOfGausspoints;
        bodyLoad.typeOfLoad = 'deadLoad';
        bodyLoad.timeFunction = @(t) t;
        bodyLoad.loadFunction = matlabFunction(formula(BSymbolic));

        %% Dirichlet BC
        dirichletSurfaceX = dirichletClass(dofObject);
        dirichletSurfaceX.nodalDof = 1;
        dirichletSurfaceX.masterObject = solidElectroThermoObject;
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
        dirichletSurfaceY.masterObject = solidElectroThermoObject;
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
        dirichletSurfaceZ.masterObject = solidElectroThermoObject;
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
        dirichletElectricalPotential.nodalDof = solidElectroThermoObject.dimension + 1;
        dirichletElectricalPotential.masterObject = solidElectroThermoObject;
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

        if choiceFunction == 2
            dirichletThermo = dirichletClass(dofObject);
            dirichletThermo.functionDof = 2;%'addfield2Y';
            dirichletThermo.nodalDof = solidElectroThermoObject.dimension + 2;
            dirichletThermo.masterObject = solidElectroThermoObject;
            dirichletThermo.nodeList = auessereKnoten;
            dirichletThermo.timeFunction = str2func('@(t,X2) 293.15 + t*10*cos(X2)');
        elseif choiceFunction == 3
            dirichletThermo = dirichletClass(dofObject);
            dirichletThermo.functionDof = 2;%'addfield2Y';
            dirichletThermo.nodalDof = solidElectroThermoObject.dimension + 2;
            dirichletThermo.masterObject = solidElectroThermoObject;
            dirichletThermo.nodeList = auessereKnoten;
            dirichletThermo.timeFunction = str2func('@(t,X2) 293.15 + t*10*X2.^3');
        elseif choiceFunction == 4
            dirichletThermo = dirichletClass(dofObject);
            dirichletThermo.functionDof = 1:3;%'addfield2XYZ';
            dirichletThermo.nodalDof = solidElectroThermoObject.dimension + 2;
            dirichletThermo.masterObject = solidElectroThermoObject;
            dirichletThermo.nodeList = auessereKnoten;
            dirichletThermo.timeFunction = str2func('@(t,X1,X2,X3) 293.15 + t*10*X1.^2 + t*10*X2.^2 + t*10*X3.^2');
        elseif choiceFunction == 5
            KnotenLinks = unique([find(nodes(:,1) == 0)]);
            KnotenRechts = unique([find(nodes(:,1) == 1)]);
            dirichletThermoLinks = dirichletClass(dofObject);
            dirichletThermoLinks.functionDof = 1;%'addfield2X';
            dirichletThermoLinks.nodalDof = solidElectroThermoObject.dimension + 2;
            dirichletThermoLinks.masterObject = solidElectroThermoObject;
            dirichletThermoLinks.nodeList = KnotenLinks;
            dirichletThermoLinks.timeFunction = str2func('@(t,X1) t*20*X1.^3');
            dirichletThermoRechts = dirichletClass(dofObject);
            dirichletThermoRechts.functionDof = 1;%'addfield2X';
            dirichletThermoRechts.nodalDof = solidElectroThermoObject.dimension + 2;
            dirichletThermoRechts.masterObject = solidElectroThermoObject;
            dirichletThermoRechts.nodeList = KnotenRechts;
            dirichletThermoRechts.timeFunction = str2func('@(t,X1) t*20*X1.^3');
        else
            error('not implemented');
        end

        %% Test
        % xanaly
        if choiceFunction == 2
            for indexX = 0:0.5:1
                if abs(double(phiFunction(indexX,0,0) - dirichletElectricalPotential.timeFunction(1,indexX))) >= 1e-12
                    error('Analytical formulations are divergent')
                end
                for indexY = 0:0.5:1
                    if abs(double(thetaFunction(0,indexY,0)-dirichletThermo.timeFunction(1,indexY))) >= 1e-12
                        error('Analytical formulations are divergent')
                    end
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
                        if abs(double(thetaFunction(indexX,indexY,indexZ)-dirichletThermo.timeFunction(1,indexX,indexY,indexZ))) >= 1e-12
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
            % temperature field
            analyticalVariables.theta = zeros(size(nodes(:,1)));
            % electrical Displacement
            analyticalVariables.D = zeros(size(nodes(:,1:3)));
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
                analyticalVariables.phi(i,:)        = phiFunction(nodes(i,1),nodes(i,2),nodes(i,3));
                analyticalVariables.theta(i,:)      = thetaFunction(nodes(i,1),nodes(i,2),nodes(i,3));
                analyticalVariables.D(i,:)          = DFunction(nodes(i,1),nodes(i,2),nodes(i,3));
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
            analyticalVariablesCell{ii,j}.q      = [analyticalVariables.x,analyticalVariables.phi,analyticalVariables.theta];
            analyticalVariablesCell{ii,j}.D      = analyticalVariables.D;
            analyticalVariablesCell{ii,j}.C      = analyticalVariables.C;
            analyticalVariablesCell{ii,j}.G      = analyticalVariables.G;
            analyticalVariablesCell{ii,j}.c      = analyticalVariables.c;
            analyticalVariablesCell{ii,j}.lambdaC = analyticalVariables.lambdaC;
            analyticalVariablesCell{ii,j}.lambdaG = analyticalVariables.lambdaG;
            analyticalVariablesCell{ii,j}.lambdac = analyticalVariables.lambdac;
            % % %
            %             thetaAnalyticalMid(j) = thetaAnalyticFunction(0.5,0.5,0.5);
            %             indexMidNode = find(sum((part.NODES==0.5),2)==3);
            %             thetaNumericalMid(j) = solidElectroThermoSystem.QN1(indexMidNode,5);
            solidElectroThermoObjectCell{ii,j} = solidElectroThermoObject;
            Elem{1}.disp(j,:) = solidElectroThermoObject.qN1(find(nodes(:,1) == 0.5 & nodes(:,2) == 0.5 & nodes(:,3) == 0.5),1:3) - 0.5;
        end
    end
end
save(strcat(setupObject.saveObject.fileName,'Cell3'), 'solidElectroThermoObjectCell', 'analyticalVariablesCell','-v7.3');
% save(strcat(setupObject.saveObject.fileName,'Cell'), 'solidElectroThermoObjectCell', 'analyticalVariablesCell', '-append');


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
