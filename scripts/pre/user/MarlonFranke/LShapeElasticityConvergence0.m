clear all;
%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'LShape';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 500;
setupObject.totalTime = 5;
setupObject.plotObject.flag = false;
setupObject.plotObject.view = [0,90];
setupObject.plotObject.colorBarLimits = [0 600];
setupObject.plotObject.border = 0.4;
setupObject.plotObject.lineWidth = 2;
setupObject.plotObject.colorscheme = 'turbo';
setupObject.plotObject.steps = 1;
setupObject.plotObject.savePlot.flag = false;
setupObject.plotObject.savePlot.name = 'LShape';
setupObject.plotObject.savePlot.flagPlotColorbar = false;
setupObject.newton.tolerance = 1e-6; % not good for ANN model<
setupObject.integrator = 'Midpoint';
% setupObject.integrator = 'DiscreteGradient';
dofObject = dofClass;   % required object for dof and object handling
%% continuum Objects
% abaqus mesh
mesh = 'H1H0';
% mesh = 'H1';
solidObject = solidClass(dofObject);
% abaqusMeshName = 'LShapeH1VeryCoarse';
abaqusMeshName = 'LShapeH1';
% abaqusMeshName = 'LShapeH1Medium';
abaqusMeshData = abaqusInputFileConverter(strcat(abaqusMeshName,'.inp'));
if strcmpi(mesh,'H1')
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    % solidObject.elementDisplacementType = 'displacementSC';
    solidObject.elementDisplacementType = 'displacementSC';
    % solidObject.elementNameAdditionalSpecification = 'pH';
elseif strcmpi(mesh,'H1H0')
    solidObject.shapeFunctionObject.order = 1;
    solidObject.shapeFunctionObject.numberOfGausspoints = 8;
    solidObject.mixedFEObject.typeShapeFunctionData = 0;
    solidObject.elementDisplacementType = 'mixedSC';
    % solidObject.elementDisplacementType = 'displacementSC';
    solidObject.flagHistoryFields = true;
    % solidObject.elementNameAdditionalSpecification = 'pHGJLambda';
    solidObject.elementNameAdditionalSpecification = 'pHCGJLambda';
    % solidObject.elementNameAdditionalSpecification = 'pHCGJ';
    solidObject.mixedFEObject.condensation = true;
end
solidObject.meshObject.nodes = abaqusMeshData.qR;
solidObject.meshObject.edof = abaqusMeshData.edof;
solidObject.numericalTangentObject.computeNumericalTangent = true;
% solidObject.numericalTangentObject.type = 'complex';
solidObject.materialObject.name = 'MooneyRivlin';     % ground truth
%
bulk    = 5209;                                                 % K
shear   = 997.5;                                                % G
mu      = shear;                                                % second lame parameter
a      = 5/6*mu;                                               % a
b      = 1/6*mu;                                               % b
c       = 10000;
solidObject.materialObject.a = a;
solidObject.materialObject.b = b;
solidObject.materialObject.c = c;
solidObject.materialObject.d = 2*(solidObject.materialObject.a + 2*solidObject.materialObject.b);
solidObject.materialObject.lambda = 10000;
solidObject.materialObject.mu = mu;
solidObject.materialObject.rho = 100;
solidObject.dimension = 3;
FA = [256;512;768];
bcTimeEnd = 5;
%
neumannObject1 = neumannClass(dofObject);
neumannObject1.loadType = 'deadLoad';
neumannObject1.loadGeometry = 'area';
neumannObject1.masterObject = solidObject;
neumannObject1.loadVector = 1/9*FA;
neumannObject1.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject1.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject1.projectionType = 'none';
neumannObject1.timeFunction = @(t) 1/(bcTimeEnd/2)*t.*(t <= bcTimeEnd/2)+1/(bcTimeEnd/2)*(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject1.meshObject.edof = abaqusMeshData.subsets(3).edof;
%
neumannObject2 = neumannClass(dofObject);
neumannObject2.loadType = 'deadLoad';
neumannObject2.loadGeometry = 'area';
neumannObject2.masterObject = solidObject;
neumannObject2.loadVector = -1/9*FA;
neumannObject2.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject2.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject2.projectionType = 'none';
neumannObject2.timeFunction =  @(t) 1/(bcTimeEnd/2)*t.*(t <= bcTimeEnd/2)+1/(bcTimeEnd/2)*(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject2.meshObject.edof = abaqusMeshData.subsets(4).edof;
%% solver
dofObject = runNewton(setupObject,dofObject);
save(strcat('lShapeLoadPeriod',solidObject.elementDisplacementType,solidObject.elementNameAdditionalSpecification,setupObject.integrator,abaqusMeshName,'.mat'))

% LShapeElasticityConvergence1