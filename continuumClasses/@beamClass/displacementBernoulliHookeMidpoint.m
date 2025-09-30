function [rData, kData, elementData, array] = displacementBernoulliHookeMidpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, ~, ~)
%BEAMBERNOULLIENDPOINT Bernoulli beam

%% SETUP
% load objects
meshObject = obj.meshObject;
materialObject = obj.materialObject;
shapeFunctionObject = obj.shapeFunctionObject;
% material data
E = materialObject.E;
I = materialObject.I;

% geometry
dimension = obj.dimension;
assert(dimension == 1, 'Beam element is defined only for 1D!');
edof = meshObject.edof;
edRef = meshObject.nodes(edof(e, :), 1:dimension).';
he = abs(edRef(1)-edRef(2)); % element length
Je = he/2;                   % Jacobian is half of element length
detJ = det(Je);     
% deflection w and slope w' at the nodes of the element
wN1 = dofs.edN1;
wApostropheN1 = dofs.phiN1;
uN1 = [wN1(1); wApostropheN1(1); wN1(2);wApostropheN1(2)];

wN = obj.qN(edof(e, :), 1);
wApostropheN = obj.qN(edof(e, :), 2);
uN = [wN(1);wApostropheN(1);wN(2);wApostropheN(2)];
uN05 = ( uN + uN1 ) / 2;
           

% gauss integration and cubic Hermite shape functions
gaussWeight = shapeFunctionObject.gaussWeight;
numberOfGausspoints = shapeFunctionObject.numberOfGausspoints;
gaussPoints = shapeFunctionObject.gaussPoint;
[~, ~, d2N_xi_xi_k_I, d3N_xi_xi_xi_k_I] = computeCubicHermiteShapeFunctions(dimension,numberOfGausspoints,gaussPoints);

% initialize residual & tangent
RW = rData{1};
RWprime = rData{2};
KWW = kData{1, 1};
KWWprime = kData{1, 2};
KWprimeW = kData{2, 1};
KWprimeWprime = kData{2, 2};

% initialize elementEnergy
elementData.strainEnergy = 0;

%stiffness matrix
stiffness_mat = zeros(4,4);

%% Gauss loop
for k = 1:numberOfGausspoints
    % B-matrix is transformed second derivatives w.r.t. physical coordinate
    B_matrix = d2N_xi_xi_k_I(k,:) .* [1, Je, 1, Je] * Je^(-2); 
    kappaN1 = B_matrix*uN1; %curvature
    % shear force computation
    B_shear_matrix = d3N_xi_xi_xi_k_I(k,:) .* [1, Je, 1, Je] * Je^(-3); 
    wTriplePrime = B_shear_matrix*uN1;

    if ~computePostData
        
        % stiffness matrix
        stiffness_mat = stiffness_mat + E*I*(B_matrix')*B_matrix * detJ * gaussWeight(k); 
        
        % residual vector (linear)
        R = stiffness_mat*uN05;

        elementData.strainEnergy = elementData.strainEnergy + 1/2 * E*I * kappaN1^2 * detJ * gaussWeight(k);
    else

        % stress resultants
        bending_moment = -E*I*kappaN1;
        shear_force = -E*I*wTriplePrime;
        stressTensor.Cauchy = [bending_moment,shear_force];
        array = postStressComputation(array,shapeFunctionObject.N_k_I,k,gaussWeight,detJ,stressTensor,setupObject,dimension);

    end
end

%% PASS COMPUTATION DATA
if ~computePostData

    % pass residual
    rData{1} = RW + [R(1);R(3)];
    rData{2} = RWprime + [R(2);R(4)];

    % pass tangent
    kData{1, 1} = KWW + 0.5*[stiffness_mat(1,1), stiffness_mat(1,3); stiffness_mat(3,1), stiffness_mat(3,3)];
    kData{1, 2} = KWWprime + 0.5*[stiffness_mat(1,2), stiffness_mat(1,4); stiffness_mat(3,2), stiffness_mat(3,4)];
    kData{2, 1} = KWprimeW + 0.5*[stiffness_mat(2,1), stiffness_mat(2,3); stiffness_mat(4,1), stiffness_mat(4,3)];
    kData{2, 2} = KWprimeWprime + 0.5*[stiffness_mat(2,2), stiffness_mat(2,4); stiffness_mat(4,2), stiffness_mat(4,4)];
end

end