function displacementViscoEndpoint(obj,setupObject,computePostData)
% homogenous linear-viscoelastic isotropic strain-energy function, evaluated at time n+1, i.e. implicid euler method.
% 08.02.2015 Marlon Franke: first esra version
% 19.08.2021 Marlon Franke: moofeKIT version
%creates the residual and the tangent of the given obj.

%% setup
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
storageFEObject = obj.storageFEObject;
materialObject = obj.materialObject;
mapVoigtObject = obj.mapVoigtObject;
%   mixedFEObject = obj.mixedFEObject;

meshObject = obj.meshObject;
% element degree of freedom tables and more
edof = meshObject.edof;
globalFullEdof = meshObject.globalFullEdof;
numberOfElements = size(globalFullEdof,1);
numberOfDOFs = size(globalFullEdof,2);
dimension = obj.dimension;

% gauss integration and shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
%   N = shapeFunctionObject.N;
dNr = shapeFunctionObject.dNr;
% nodal dofs
qR = obj.qR;
qN1 = obj.qN1;

% material data and voigt notation
lambda = materialObject.lambda;
mu = materialObject.mu;
eModul0 = materialObject.eModul0;
eModul1 = materialObject.eMoodul1;
eta1 = materialObject.eta1;

DT = setupObject.timeStepSize;

selectMapVoigt(mapVoigtObject,dimension,'symmetric');

if (dimension == 1)
    DMat=2*mu;
else
    DMat = [lambda+2*mu lambda lambda 0 0 0;...
        lambda lambda+2*mu lambda 0 0 0;...
        lambda lambda lambda+2*mu 0 0 0;...
        0 0 0 mu 0 0;...
        0 0 0 0 mu 0;...
        0 0 0 0 0 mu];
end

% initialization energy and storageFEObject
timeStep = setupObject.timeStep;
obj.elementData(timeStep).strainEnergy = 0;
initializeDataFE(storageFEObject);
%% element loop for residual and tangent
for e = 1:numberOfElements
    array = struct('Re', zeros(numberOfDOFs,1), 'Ke', zeros(numberOfDOFs),...% solver: computing residual and tangent
        'Se',zeros(numberOfDOFs/dimension,1),'Me',zeros(numberOfDOFs/dimension)); % postprocessing: computing stress etc.
    %post processing stresses
    edN1 = qN1(edof(e,:),1:dimension)';
    edRef = qR(edof(e,:),1:dimension)';
    uN1 = edN1 - edRef;
    J = qR(edof(e,:),1:dimension)'*dNr';

    %C = eModul0 + eModul1/(1+DT/(eta1/eModul1));
    %% Gauss loop
    for k = 1:numberOfGausspoints
        indx = dimension*k-(dimension-1):dimension*k;
        detJ = det(J(:,indx)');
        if detJ < 10*eps
            error('Jacobi determinant equal or less than zero.')
        end
        dNx = (J(:,indx)')\dNr(indx,:);
        B = BMatrix(dNx,'mapVoigtObject',mapVoigtObject);

        if ~computePostData
            %epsilonGausspoint = epsilonGaussponit + B*uN1(:);
            %sigmaGausspoint = sigmaGausspoint + (eModul0+eModul1/(1+DT/(eta1/eModul1)))*epsilonGausspoint-eModul1/(1+DT/(eta1/eModul1))*epsilon;

            % array.Re = array.Re + B'*sigmaGausspoint*detJ*gaussWeight(k);
            % array.Ke = array.Ke + B'*C*B*detJ*gaussWeight(k);

            % Residual
            array.Re = array.Re + (B'*DMat*B)*uN1(:)*detJ*gaussWeight(k);
            % Tangent
            array.Ke = array.Ke + B'*DMat*B*detJ*gaussWeight(k);
            %strain energy
            obj.elementData(timeStep).strainEnergy = obj.elementData(timeStep).strainEnergy + 1/2*(B*uN1(:))'*DMat*(B*uN1(:))*detJ*gaussWeight(k);
        else
            % stress at gausspoint
            SN1_v = DMat*EN1_v;
            SN1 = voigtToMatrix(SN1_v, 'stress');
            PN1 = FN1*SN1;
            stressTensor.FirstPK = PN1;
            stressTensor.Cauchy = 1/det(FN1)*PN1*FN1';
            array = postStressComputation(array,N,k,gaussWeight,detJ,detJN1,stressTensor,setupObject,dimension);
        end
    end
    storeDataFE(storageFEObject,obj,array,globalFullEdof,e,dimension,computePostData);
end
end