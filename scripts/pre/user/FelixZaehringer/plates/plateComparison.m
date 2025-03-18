displacementElementFormulationMatrix=zeros(7,max(size(0.0:0.1:0.4)));
formulationName=strings(8,1);
formulationName(1,:)="pure displacement 4 node displacement element";
formulationName(2,:)='displacement 4 node, SRI element';
formulationName(3,:)='displacement 4 node Bathe/Dvorkin element';
formulationName(4,:)='EAS 4 node Andelfinger-Ramm element';
formulationName(5,:)='EAS 4 node Simo-Rifai element';
formulationName(6,:)='displacement 4 node Petrov Galerkin element';
formulationName(7,:)='EAS Andelfinder-Ramm 4 node Petrov Galerkin element';
formulationName(8,:)='EAS Simo-Rifai 4 node Petrov Galerkin element';
formulationName(9,:)='EAS Simo-Rifai 4 node Petrov Galerkin element ANS';


for counter=1:3
    counter2=0;
    for ii=0.0:0.1:0.4
        counter2=counter2+1;
%% setup (mandatory: setup and dofs)
setupObject = setupClass;
%setupObject.saveObject.fileName = 'plateSimple';
setupObject.saveObject.saveData = false;
setupObject.totalTimeSteps = 1;
setupObject.totalTime = 1;
setupObject.plotObject.flag = false;
setupObject.plotObject.postPlotType = 'zero';
setupObject.newton.tolerance = 1e-6;

dofObject = dofClass;   % required object for dof and object handling

%% continuum Objects
plateLengthX = 6;
plateLengthY = 9;
orderShapeFunctions =1;
numberOfElementsX=10;
% if counter==6
% numberOfElementsX=100;
% end
numberOfElementsY=numberOfElementsX;
plateObject = plateClass(dofObject);
[plateObject.meshObject.nodes,plateObject.meshObject.edof] = meshRectangle(plateLengthX, plateLengthY, numberOfElementsX, numberOfElementsY, orderShapeFunctions);
hX=plateLengthX/(numberOfElementsX);
hY=plateLengthY/(numberOfElementsY);
h=max(hX,hY);
deplacedNode1=numberOfElementsX.*numberOfElementsY/2;
deplacedNode2=1*numberOfElementsX.*numberOfElementsY/4;
deplacedNode3=3*numberOfElementsX.*numberOfElementsY/4;
plateObject.meshObject.nodes(deplacedNode1,1)=plateObject.meshObject.nodes(deplacedNode1,1)-h*ii;
plateObject.meshObject.nodes(deplacedNode1,2)=plateObject.meshObject.nodes(deplacedNode1,2)+h*ii;
plateObject.meshObject.nodes(deplacedNode2,1)=plateObject.meshObject.nodes(deplacedNode2,1)-h*ii;
plateObject.meshObject.nodes(deplacedNode2,2)=plateObject.meshObject.nodes(deplacedNode2,2)+h*ii;
plateObject.meshObject.nodes(deplacedNode3,1)=plateObject.meshObject.nodes(deplacedNode3,1)-h*ii;
plateObject.meshObject.nodes(deplacedNode3,2)=plateObject.meshObject.nodes(deplacedNode3,2)-h*ii;





mesh=zeros(max(size(plateObject.meshObject.nodes)),2);
mesh(:,1)=plateLengthX/2;
mesh(:,2)=plateLengthY/2;
plateObject.meshObject.nodes = plateObject.meshObject.nodes + mesh;
[plateObject.elementDisplacementType, plateObject.mixedFEObject.condensation,plateObject.elementNameAdditionalSpecification,plateObject.materialObject.name,plateObject.mixedFEObject.typeShapeFunctionData] = getFormulation(counter);
plateObject.materialObject.rho = 0;
plateObject.materialObject.E = 1000;
plateObject.materialObject.nu = 0.2;
plateObject.h = 0.2;
plateObject.dimension = 2;
plateObject.shapeFunctionObject.order = orderShapeFunctions;
plateObject.shapeFunctionObject.numberOfGausspoints = (orderShapeFunctions+1)^2;




disp(formulationName(counter,:))


% Dirichlet boundary
boundary1 = dirichletClass(dofObject);
boundary1.nodeList = find(plateObject.meshObject.nodes(:,1)== 0| plateObject.meshObject.nodes(:,1) == plateLengthX | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLengthY);
boundary1.nodalDof = 1;
boundary1.masterObject = plateObject;


boundary2 = dirichletClass(dofObject);
boundary2.nodeList = find(plateObject.meshObject.nodes(:,1)== 0| plateObject.meshObject.nodes(:,1) == plateLengthX | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLengthY);
boundary2.nodalDof = 2;
boundary2.masterObject = plateObject;

boundary3 = dirichletClass(dofObject);
boundary3.nodeList = find(plateObject.meshObject.nodes(:,1)== 0| plateObject.meshObject.nodes(:,1) == plateLengthX | plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLengthY);
boundary3.nodalDof = 3;
boundary3.masterObject = plateObject;



% Neumann boundary
neumannObject = neumannClass(dofObject);
neumannObject.dimension = 3;                                % Dimension des zu betrachtenden Körpers angeben
neumannObject.masterObject = plateObject;
neumannObject.forceVector = -[100 ; 500; 200];                         % Komponenten des Kraftvektors anpassen
neumannObject.shapeFunctionObject.dimension = plateObject.dimension;
neumannObject.shapeFunctionObject.order = plateObject.shapeFunctionObject.order;
neumannObject.shapeFunctionObject.numberOfGausspoints = plateObject.shapeFunctionObject.numberOfGausspoints;  % Anzahl der Gaußpunkten auf dem Neumannrand anpassen
neumannObject.meshObject.edof = plateObject.meshObject.edof;

% % Nodal Forces
% nodalForceObject = nodalForceClass(dofObject);
% nodalForceObject.dimension = 3;                                % Dimension des zu betrachtenden Körpers angeben
% nodalForceObject.masterObject = plateObject;
% nodalForceObject.forceVector = -[0; 0.5; 0];                         % Komponenten des Kraftvektors anpassen
% nodalForceObject.nodeList = find(plateObject.meshObject.nodes(:,1) == plateLength & (plateObject.meshObject.nodes(:,2) == 0 | plateObject.meshObject.nodes(:,2) == plateLength));

%% solver
dofObject = runNewton(setupObject,dofObject);
%plot(plateObject, setupObject)
disp(['Maximum vertical displacement: ', num2str(max(abs(plateObject.qN1(:, 1))))]);
displacementElementFormulationMatrix(counter,counter2)=max(abs(plateObject.qN1(:, 1)));
    end
end
%% Postprocessing


MeshDistortion=figure;
displacementElementFormulationMatrix2=displacementElementFormulationMatrix;
x=0.0:0.1:0.4;
for counter=1:3
    plot(x,displacementElementFormulationMatrix(counter,:), 'DisplayName', formulationName(counter,:))
    displacementElementFormulationMatrix2(counter,:)=100*displacementElementFormulationMatrix2(counter,:)/displacementElementFormulationMatrix(counter,1);
hold on
hl = legend('show');
ylabel('displacement $w$','Interpreter','Latex')
xlabel('mesh distortion $s$','Interpreter','Latex')
set(hl, 'Interpreter','latex')
grid()
title('Mesh distortion test')
%exportgraphics(MeshDistortion,'MeshDistortion.eps')
%exportgraphics(MeshDistortion,'MeshDistortion.jpg')

end

hold off

MeshDistortion2=figure;
x=0.0:0.1:0.4;
for counter=1:3
    plot(x,displacementElementFormulationMatrix2(counter,:), 'DisplayName', formulationName(counter,:))
hold on
hl = legend('show');
ylabel('relative displacement $\%$','Interpreter','Latex')
xlabel('mesh distortion $s$','Interpreter','Latex')
set(hl, 'Interpreter','latex')
grid()
title('Mesh distortion test')
%exportgraphics(MeshDistortion2,'MeshDistortion2.eps')
%exportgraphics(MeshDistortion2,'MeshDistortion2.jpg')

end



function [elementDisplacementType, condensation,elementNameAdditionalSpecification,materialName,typeShapeFunctionData] = getFormulation(counter)

typeShapeFunctionData='';
% 1. pure displacement 4 node displacement element
% 2. displacement 4 node, SRI element
% 3. displacement 4 node Bathe/Dvorkin element
% 4. EAS Andelfinger-Ramm element
% 5. EAS Simo-Rifai element
% 6. displacement 4 node Petrov Galerkin element
% 7. EAS Andelfinder-Ramm 4 node Petrov Galerkin element
% 8. EAS Simo-Rifai 4 node Petrov Galerkin element 
% 9. EAS Simo-Rifai 4 node Petrov Galerkin element Ans 

condensation = false;
elementNameAdditionalSpecification='';
if counter==1
    elementDisplacementType='displacement';
    materialName='Hooke';

elseif counter==2
%     elementDisplacementType='selectiveReducedIntegration';
    materialName='Hooke';
    elementDisplacementType='displacement';
%     condensation = true;

elseif counter==3
    elementDisplacementType='displacement';
    materialName='HookeBD';


elseif counter==4
    elementDisplacementType='eas';
    materialName='HookeBD';
%     elementNameAdditionalSpecification='AndelfingerRamm';
    typeShapeFunctionData = 4;


elseif counter==5
    elementDisplacementType='eas';
    materialName='Hooke';
    elementNameAdditionalSpecification='SimoRifai';
   typeShapeFunctionData = 4;


elseif counter==6
    elementDisplacementType = 'displacementPetrovGalerkin';
    materialName='Hooke';

elseif counter==7
    elementDisplacementType = 'easPetrovGalerkin';
    materialName='HookeBD';
    typeShapeFunctionData = 4;
%     elementNameAdditionalSpecification='AndelfingerRamm';
elseif counter==8
    elementDisplacementType = 'easPetrovGalerkin';
    materialName='HookeBD';
    typeShapeFunctionData = 4;
%     elementNameAdditionalSpecification='SimoRifai';
elseif counter==9
    elementDisplacementType = 'easPetrovGalerkin';
    materialName='HookeBD';
    typeShapeFunctionData = 4;
%     elementNameAdditionalSpecification='SimoRifai';




end

end
