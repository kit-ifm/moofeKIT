% LSHAPE Script for preprocessing a dynamic mechanical simulation.
% 
% FORMULATION
% Different formulations like standard displacement-based and mixed en-
% hanced assumed strain (eas) and different material models can be chosen.
% 
% REFERENCE
% https://doi.org/10.1007/BF00913408
% 
% SEE ALSO 
% cooksMembrane,
% LShapeElectroThermo
% 
% CREATOR(S) 
% Marlon Franke

%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'LShape';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 80;
setupObject.totalTime = 8;
setupObject.plotObject.flag = false;
% setupObject.plotObject.stress.component = -1;
setupObject.plotObject.view = [-0.4,90];
setupObject.newton.tolerance = 1e-4;
setupObject.integrator = 'Endpoint';
% setupObject.integrator = 'Midpoint';
% setupObject.integrator = 'DiscreteGradient';

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
% abaqus mesh
% abaqusMeshData = abaqusInputFileConverter('LShapeH1.inp');
abaqusMeshData = abaqusInputFileConverter('LShapeH2Serendipity.inp');
solidObject = solidClass(dofObject);
solidObject.meshObject.nodes = abaqusMeshData.qR;
solidObject.meshObject.edof = abaqusMeshData.edof;
if 1
    solidObject.materialObject.name = 'MooneyRivlin';
%     solidObject.elementDisplacementType = 'displacementSC';
    solidObject.elementDisplacementType = 'mixedSC';
    bulk    = 5209;                                                 % K
    shear   = 997.5;                                                % G
    mu      = shear;                                                % second lame parameter
    c1      = 5/6*mu;                                               % a
    c2      = 1/6*mu;                                               % b
    c       = 0;
    solidObject.materialObject.a = c1/2;
    solidObject.materialObject.b = c2/2;
    solidObject.materialObject.c = c;
    solidObject.materialObject.d = 2*(solidObject.materialObject.a + 2*solidObject.materialObject.b);
    solidObject.mixedFEObject.condensation = true;
    solidObject.mixedFEObject.typeShapeFunctionData = 1;
elseif 0
    solidObject.materialObject.name = 'NeoHooke';
%     solidObject.materialObject.name = 'SaintVenant';
%     solidObject.materialObject.name = 'SaintVenantNumericalTangent';
    % solidObject.materialObject.name = 'Hyperelastic'; % 'Hooke';
    solidObject.elementDisplacementType = 'displacement'; %'displacement'; % 'displacementSC'; % 'eas';
%     solidObject.elementDisplacementType = 'displacementSC';
%     solidObject.elementDisplacementType = 'eas';
    solidObject.mixedFEObject.condensation = true;
%     solidObject.mixedFEObject.condensation = true;
    % 'displacementSC'; % 'eas';
    solidObject.materialObject.rho = 1;
    E = 2.26115e+03;
    nu = 4.95469e-01;
    solidObject.materialObject.lambda = E*nu/((1+nu)*(1-2*nu));
    solidObject.materialObject.mu = E/(2*(1+nu));                              % G (second lame parameter)
elseif 0
    solidObject.materialObject = materialPackage.saintVenantClass;
%     solidObject.materialObject.name = 'Hyperelastic'; % 'Hooke';
%     solidObject.elementDisplacementType = 'displacement'; % 'displacementSC'; % 'eas';
    solidObject.materialObject.name = 'Hooke';
    solidObject.elementDisplacementType = 'displacement'; % 'displacementSC'; % 'eas';
%     solidObject.elementDisplacementType = 'eas';
    solidObject.materialObject.E = 2.26115e+03;
    solidObject.materialObject.nu = 4.95469e-01;
elseif 0
    solidObject.materialObject = materialPackage.mooneyRivlinClass;
    solidObject.materialObject.name = 'Hyperelastic'; % 'Hooke';
    solidObject.elementDisplacementType = 'displacement'; % 'displacementSC'; % 'eas';
    bulk    = 5209;                                                 % K
    shear   = 997.5;                                                % G
    mu      = shear;
    solidObject.materialObject.c1 = 5/6*mu;
    solidObject.materialObject.c2 = 1/6*mu;
    solidObject.materialObject.c = 0;    
end
solidObject.materialObject.rho = 10;
solidObject.dimension = 3;
solidObject.shapeFunctionObject.order = 2;
% solidObject.shapeFunctionObject.numberOfGausspoints = 8;
solidObject.shapeFunctionObject.numberOfGausspoints = 27;
% solidObject.numericalTangentObject.computeNumericalTangent = true;
% solidObject.numericalTangentObject.showDifferences = true;

bcTimeEnd = 5;
FA = [256;512;768];
%
neumannObject1 = neumannClass(dofObject);
neumannObject1.typeOfLoad = 'deadLoad';
neumannObject1.masterObject = solidObject;
neumannObject1.forceVector = 1/9*FA;
neumannObject1.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject1.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject1.projectionType = 'none';
neumannObject1.timeFunction = @(t) t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject1.meshObject.edof = abaqusMeshData.subsets(3).edof;
%
neumannObject2 = neumannClass(dofObject);
neumannObject2.typeOfLoad = 'deadLoad';
neumannObject2.masterObject = solidObject;
neumannObject2.forceVector = -1/9*FA;
neumannObject2.shapeFunctionObject.order = solidObject.shapeFunctionObject.order;
neumannObject2.shapeFunctionObject.numberOfGausspoints = 2^(solidObject.dimension-1);
neumannObject2.projectionType = 'none';
neumannObject2.timeFunction = @(t)  t.*(t <= bcTimeEnd/2)+(bcTimeEnd-t).*(t > bcTimeEnd/2 & t <= bcTimeEnd);
neumannObject2.meshObject.edof = abaqusMeshData.subsets(4).edof;

%% solver
dofObject = runNewton(setupObject,dofObject);

%% postprocessing - energy
timeVector = getTime(dofObject.postDataObject,setupObject);
kineticEnergy = getKineticEnergy(dofObject.postDataObject,setupObject);
[linearMomentum, totalLinearMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'L',3);
[angularMomentum, totalAngularMomentum] = getMomentum(dofObject.postDataObject,dofObject,setupObject,'J',3);
strainEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'strainEnergy');
externalEnergy = getEnergy(dofObject.postDataObject,dofObject,setupObject,'externalEnergy');
totalEnergy = strainEnergy + kineticEnergy;
tStartDiff = ceil(bcTimeEnd/setupObject.totalTime*setupObject.totalTimeSteps);
totalEnergyDiff = totalEnergy(tStartDiff+1:end) - totalEnergy(tStartDiff:end-1);
figure;
plot(timeVector, totalEnergy);
figure;
plot(timeVector(tStartDiff+1:end), totalEnergyDiff);
matlab2tikz(['diagram','.tikz'],'width','\figW','height','\figH')
figure;
plot(timeVector,linearMomentum);
figure;
plot(timeVector,totalLinearMomentum);
figure;
plot(timeVector,angularMomentum);
figure;
plot(timeVector,totalAngularMomentum);
