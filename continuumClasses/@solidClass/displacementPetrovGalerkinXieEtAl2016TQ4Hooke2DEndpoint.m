function [rData, kData, elementEnergy, array] = displacementPetrovGalerkinXieEtAl2016TQ4Hooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% DISPLACEMENTPETROVGALERKINXIEETAL2016TQ4HOOKE2DENDPOINT Element routine of class solidClass.
%
% FORMULATION
% This is a 'displacement'-based finite element routine covering linear
% mechanical processes employing a homogenous, linear-elastic, isotropic
% 'Hooke' material model (linear geometric and linear material/
% stress-strain relation).
% The formulation employs different shape functions for trial and test
% function ('PetrovGalerkin') and is suitable for static simulations.
% Another name of this element is TQ4, where T denotes the Trefftz unsymmetry
% of the stiffness matrix.
%
% CALL
% displacementPetrovGalerkinXieEtAl2016TQ4Hooke2DEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, stressTensor, flagNumericalTangent)
% obj: The first argument is expected to be an object of type solidClass,
%      e.g. solidObject.
% setupObject: The second argument is expected to be an object of type
%              setupClass, e.g. setupObject which contains information like
%              time step size or plotting information.
% computePostData: Logical data type which is true for computing stress
%                  only and false for computing residual and tangent.
% e: current element number
% rData: cell-array of size [1, 1] for residual data.
% kData: cell-array of size [1, 1] for tangent data.
% dofs: degrees of freedom (dofs) optionally manipulated data (numerical
%       tangent)
% array: structure for storage fe-data, for more information see
%        storageFEObject.initializeArrayStress
% stressTensor: structure for storage stress tensors (postprocessing), for
%               more information see storageFEObject.initializeArrayStress
% flagNumericalTangent: flag that indicates whether the function call
%                       happens during the computation of the numerical
%                       tangent or not.
%
% REFERENCE
% https://doi.org/10.1007/s10999-014-9289-3 (Xie, Sze, Zhou: Modified and Trefftz unsymmetric finite element models)
%
% SEE ALSO
% -
%
% CREATOR(S)
% Jakob Hammes, Felix Zaehringer

%% SETUP
% load objects
shapeFunctionObject = obj.shapeFunctionObject;
materialObject = obj.materialObject;
meshObject = obj.meshObject;

% aquire shape functions
N_k_I = shapeFunctionObject.N_k_I;

% acquire the derivatives of the shape functions with respect to the parametric coordinates
dN_xi_k_I = shapeFunctionObject.dN_xi_k_I;
dN0_xi_I = shapeFunctionObject.dN0_xi_I;

% aquire the edof (element degree of freedom table)
edof = meshObject.edof;

% aquire the number of dimensions of the simulation
dimension = obj.dimension;

% number of Gausspoints
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;

% gauss weights
gaussWeight = shapeFunctionObject.gaussWeight;
gaussPoints = shapeFunctionObject.gaussPoint;

% aquire material data
nu = materialObject.nu;
E = materialObject.E;

% compute material matrix
switch lower(strtok(obj.materialObject.name, 'Hooke'))
    case 'esz'
        % ESZ
        C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
        EBar = E;           % (necessary) for Trefftz-trial-functions                                                
        nuBar = nu;         % necessary for Trefftz-trial-functions
    case 'evz'
        % EVZ
        C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
        EBar = E/(1-nu);    
        nuBar = nu/(1-nu);  
    otherwise
        error('not implemented')
end

% spatial coordinates of the element
x = dofs.edN1;
% material coordinates of the element
X = obj.qR(edof(e, :), 1:dimension)';

% compute Jacobian matrices
JAll = computeJacobianForAllGausspoints(X, dN_xi_k_I);
J0 = X * dN0_xi_I';

% displacement vector
uN1 = x(:) - X(:);

% initialize residual
Re = rData{1};

% initialize tangent
Ke = kData{1, 1};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;


%% Computing Trefftz-trial-functions
% parameters for isoparametric coordinate transformation
parameterOfIsoparametricCoordinateTransformation = 0.25*[ 1  1  1  1;      
                                                          -1  1  1 -1;
                                                           1 -1  1 -1;
                                                          -1 -1  1  1]*X';

% name change for better readability 
aXi  = parameterOfIsoparametricCoordinateTransformation(2,1);     
aEta = parameterOfIsoparametricCoordinateTransformation(4,1);    
bXi  = parameterOfIsoparametricCoordinateTransformation(2,2);
bEta = parameterOfIsoparametricCoordinateTransformation(4,2);

% transformation matrices
tXi  = 1/(sqrt(aXi^2  + bXi^2))*[aXi -bXi;  bXi aXi];
tEta = 1/(sqrt(aEta^2 + bEta^2))*[aEta -bEta;  bEta aEta];

% nodes and Gauss-Points in x-Hat coordinates
xHat           = X-parameterOfIsoparametricCoordinateTransformation(1,:)';                      % xHat  = [xHat;yHat] = [x-x0;y-y0]
xHatGaussPoint = (parameterOfIsoparametricCoordinateTransformation(2,:)'*gaussPoints(1,:) ...   % xHatGaussPoint = [xHatGaussPoint;yHatGaussPoint]
                + parameterOfIsoparametricCoordinateTransformation(4,:)'*gaussPoints(2,:) ...
                + parameterOfIsoparametricCoordinateTransformation(3,:)'*gaussPoints(1,:).*gaussPoints(2,:));    

% nodes and Gauss-Points in rXi-sXi / rEta-sEta coordinates 
rXi  = 1/sqrt(aXi^2 + bXi^2)*[aXi bXi;-bXi aXi]*xHat;                   % rXi  = [rXi;sXi]
rEta = 1/sqrt(aEta^2 + bEta^2)*[aEta bEta;-bEta aEta]*xHat;             % rEta = [rEta;sEta]

rXiGaussPoint  = 1/sqrt(aXi^2 + bXi^2)*[aXi bXi;-bXi aXi]*xHatGaussPoint;
rEtaGaussPoint = 1/sqrt(aEta^2 + bEta^2)*[aEta bEta;-bEta aEta]*xHatGaussPoint;


% computing filled pQ-matrix (right side)
pQAllNodes = [];
% loop over all corner-nodes                
for k = 1: 4
    pQNode = [1 0 xHat(1,k) 0 xHat(2,k) 0 tXi(1,:)*[(-2*rXi(1,k)*rXi(2,k)) ; (nuBar*rXi(2,k)^2 +rXi(1,k)^2)] tEta(1,:)*[-2*rEta(1,k)*rEta(2,k) ; (nuBar*rEta(2,k)^2 + rEta(1,k)^2)] 
              0 1 0 xHat(1,k) 0 xHat(2,k) tXi(2,:)*[(-2*rXi(1,k)*rXi(2,k)) ; (nuBar*rXi(2,k)^2 +rXi(1,k)^2)] tEta(2,:)*[-2*rEta(1,k)*rEta(2,k) ; (nuBar*rEta(2,k)^2 + rEta(1,k)^2)] ];
    pQAllNodes=vertcat(pQAllNodes,pQNode);
end

% derivatives
pQdx = zeros(2,8,4);
pQdy = zeros(2,8,4);
pQGaussPointdx = zeros(2,8,4);
pQGaussPointdy = zeros(2,8,4);

% derivatives of rXi-sXi and rEta-sEta coordinates are constant
rXidx  = 1/(sqrt(aXi^2 +bXi^2))*[ aXi ; -bXi];
rXidy  = 1/(sqrt(aXi^2 +bXi^2))*[ bXi ;  aXi];
rEtadx = 1/(sqrt(aEta^2 +bEta^2))*[ aEta ; -bEta];
rEtady = 1/(sqrt(aEta^2 +bEta^2))*[ bEta ;  aEta];


% derivatives of pQ for each gauss point and computing derivatives of the Trefftz-trial-functions
% loop over each gausspoint
for nn = 1:4
    %derivatives of pQ
    pQdx(:,:,nn) = [0 0 1 0 0 0 tXi(1,:)*[-2*rXidx(1)*rXiGaussPoint(2,nn)-2*rXidx(2)*rXiGaussPoint(1,nn); 2*nuBar*rXidx(2)*rXiGaussPoint(2,nn)+2*rXidx(1)*rXiGaussPoint(1,nn)] tEta(1,:)*[-2*rEtadx(1)*rEtaGaussPoint(2,nn)-2*rEtadx(2)*rEtaGaussPoint(1,nn); 2*nuBar*rEtadx(2)*rEtaGaussPoint(2,nn)+2*rEtadx(1)*rEtaGaussPoint(1,nn)]
                    0 0 0 1 0 0 tXi(2,:)*[-2*rXidx(1)*rXiGaussPoint(2,nn)-2*rXidx(2)*rXiGaussPoint(1,nn); 2*nuBar*rXidx(2)*rXiGaussPoint(2,nn)+2*rXidx(1)*rXiGaussPoint(1,nn)] tEta(2,:)*[-2*rEtadx(1)*rEtaGaussPoint(2,nn)-2*rEtadx(2)*rEtaGaussPoint(1,nn); 2*nuBar*rEtadx(2)*rEtaGaussPoint(2,nn)+2*rEtadx(1)*rEtaGaussPoint(1,nn)]];

    pQdy(:,:,nn) = [0 0 0 0 1 0 tXi(1,:)*[-2*rXidy(1)*rXiGaussPoint(2,nn)-2*rXidy(2)*rXiGaussPoint(1,nn); 2*nuBar*rXidy(2)*rXiGaussPoint(2,nn)+2*rXidy(1)*rXiGaussPoint(1,nn)] tEta(1,:)*[-2*rEtady(1)*rEtaGaussPoint(2,nn)-2*rEtady(2)*rEtaGaussPoint(1,nn); 2*nuBar*rEtady(2)*rEtaGaussPoint(2,nn)+2*rEtady(1)*rEtaGaussPoint(1,nn)]
                    0 0 0 0 0 1 tXi(2,:)*[-2*rXidy(1)*rXiGaussPoint(2,nn)-2*rXidy(2)*rXiGaussPoint(1,nn); 2*nuBar*rXidy(2)*rXiGaussPoint(2,nn)+2*rXidy(1)*rXiGaussPoint(1,nn)] tEta(2,:)*[-2*rEtady(1)*rEtaGaussPoint(2,nn)-2*rEtady(2)*rEtaGaussPoint(1,nn); 2*nuBar*rEtady(2)*rEtaGaussPoint(2,nn)+2*rEtady(1)*rEtaGaussPoint(1,nn)]];
     
    % derivatives of Trefftz-trial-functions
    pQGaussPointdx(:,:,nn)=pQdx(:,:,nn)/pQAllNodes;
    pQGaussPointdy(:,:,nn)=pQdy(:,:,nn)/pQAllNodes;

end

% assembling T-matrices
TMatrix = [pQGaussPointdx(1,:,:); pQGaussPointdy(2,:,:) ; pQGaussPointdy(1,:,:)+pQGaussPointdx(2,:,:)];

%% GAUSS LOOP
for k = 1:numberOfGausspoints
    % compute the Jacobian determinant
    [J, detJ] = extractJacobianForGausspoint(JAll, k, setupObject, dimension);
    dN_X_k_I = computedN_X_I(dN_xi_k_I, J, k);

    % compute nodal operator matrix
    B = BMatrix(dN_X_k_I);    
    
    % extracting 'B-Matrix' of Trefftz-trial-functions
    T = TMatrix(:,:,k);

    if ~computePostData
        % TANGENT
        Ke = Ke + B' * C * T * detJ * gaussWeight(k);

        % ENERGY
        elementEnergy.strainEnergy = elementEnergy.strainEnergy + 1 / 2 * (T * uN1)' * C * (T * uN1) * detJ * gaussWeight(k);
    else
        % STRESS COMPUTATION
        sigmaN1_v = C * T * uN1;
        stressTensor.Cauchy = voigtToMatrix(sigmaN1_v, 'stress');
        array = postStressComputation(array, N_k_I, k, gaussWeight, detJ, stressTensor, setupObject, dimension);
    end
end

%% RESIDUAL
Re = Re + Ke * uN1;

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = Re;

    % pass tangent
    kData{1, 1} = Ke;
end
end