%% setup (mandatory: setup and dofs)
setupObject = setupClass;
setupObject.saveObject.fileName = 'patchTest';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 10;
setupObject.totalTime = 1;
setupObject.plotObject.flag = true;
setupObject.plotObject.stress.component = -1;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
solidObject = solidClass(dofObject);
[solidObject.meshObject.nodes,solidObject.meshObject.edof] = meshGeneratorCube(1,1,1,2,2,2,2,true);
solidObject.meshObject.nodes = solidObject.meshObject.nodes + 0.5;
% solidObject.materialObject.name = 'Hooke'; 
solidObject.materialObject.name = 'SaintVenant';
solidObject.elementDisplacementType = 'displacementSC';
solidObject.materialObject.rho = 0;
solidObject.materialObject.lambda = 27.7777;
solidObject.materialObject.mu = 41.6666;
solidObject.dimension = 3;
solidObject.shapeFunctionObject.order = 2;
solidObject.shapeFunctionObject.numberOfGausspoints = 27;

boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(solidObject.meshObject.nodes(:,3) == 0);
boundary1.nodalDof = 3;
boundary1.masterObject = solidObject;

boundary1(2) = dirichletClass(dofObject);
boundary1(2).nodeList = find(solidObject.meshObject.nodes(:,1) == 0);
boundary1(2).nodalDof = 1;
boundary1(2).masterObject = solidObject;

boundary1(3) = dirichletClass(dofObject);
boundary1(3).nodeList = find(solidObject.meshObject.nodes(:,2) == 0);
boundary1(3).nodalDof = 2;
boundary1(3).masterObject = solidObject;

boundary2 = dirichletClass(dofObject);
boundary2.nodeList = find(solidObject.meshObject.nodes(:,3) == 1);
boundary2.nodalDof = 3;
boundary2.masterObject = solidObject;
boundary2.timeFunction = str2func('@(t,Z) (Z - 1).*(t >= 1) + (Z - 1*t).*(t >= 0).*(t < 1)');

%% solver
dofObject = runNewton(setupObject,dofObject);

%% post
plot(solidObject,setupObject)
% print(figPlot,filename,'-depsc')
% savefig(figPlot,filename)
export_fig(figPlot,fileName,'-eps')
% additional code for standalone
% extraCode = {'\newcommand{\picfontsize}{\small}';
%     '\newlength\figH';
%     '\newlength\figW';
%     '\setlength{\figW}{10cm}';
%     '\setlength{\figH}{0.75\figW}'};
% tikz3dPlot(figPlot,fileName,'eps','-painters','width','\figW','floatformat','%.4g','standalone',true,...
%     'extraTikzpictureOptions','font=\picfontsize','extracode',extraCode,'showInfo',false,...
%     'extraPlotCode','\draw[blue] (0,0,0) -- (10,10,10);',...
%     'extraPlotCodeEnd',{'\draw[red] (10,0,0) circle (2pt);','\draw[red] (0,3.013,9.273) circle (2pt);',...
%     '\draw[red] (3.013,0,9.273) circle (2pt);'})


