function [rData, kData, elementEnergy, array] = beamTimoshenkoEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, ~, ~)
%BEAMTIMOSHENKOENDPOINT Timoshenko beam

%% SETUP
% load objects
meshObject = obj.meshObject;
materialObject = obj.materialObject;

% w and phi at the nodes of the element
wN1 = dofs.edN1;
phiN1 = dofs.phiN1;

% material data
E = materialObject.E;
G = materialObject.G;
I = materialObject.I;
A = materialObject.A;

% element length
dimension = obj.dimension;
assert(dimension == 1, 'Beam elements are defined only for 1D!');
edof = meshObject.edof;
nodes = meshObject.nodes(edof(e, :), 1:dimension).';
he = abs(nodes(1)-nodes(2));

% initialize residual & tangent
RW = rData{1};
RPhi = rData{2};
KWW = kData{1, 1};
KWPhi = kData{1, 2};
KPhiW = kData{2, 1};
KPhiPhi = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% TANGENT
% TODO Jonas Jundt

%% RESIDUAL
RW = RW + KWW * wN1(:) + KWPhi * phiN1(:);
RPhi = RPhi + KPhiW * wN1(:) + KPhiPhi * phiN1(:);

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = RW;
    rData{2} = RPhi;

    % pass tangent
    kData{1, 1} = KWW;
    kData{1, 2} = KWPhi;
    kData{2, 1} = KPhiW;
    kData{2, 2} = KPhiPhi;
end
end
