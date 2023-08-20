function [rData, kData, elementEnergy, array] = beamBernoulliEndpoint(obj, setupObject, computePostData, e, rData, kData, dofs, array, ~, ~)
%BEAMBERNOULLIENDPOINT Bernoulli beam

%% SETUP
% load objects
meshObject = obj.meshObject;
materialObject = obj.materialObject;

% w and w' at the nodes of the element
wN1 = dofs.edN1;
wApostropheN1 = dofs.phiN1;

% material data
E = materialObject.E;
I = materialObject.I;

% element length
dimension = obj.dimension;
assert(dimension == 1, 'Beam elements are defined only for 1D!');
edof = meshObject.edof;
nodes = meshObject.nodes(edof(e, :), 1:dimension).';
he = abs(nodes(1)-nodes(2));

% initialize residual & tangent
RW = rData{1};
RWApostrophe = rData{2};
KWW = kData{1, 1};
KWWApostrophe = kData{1, 2};
KWApostropheW = kData{2, 1};
KWApostropheWApostrophe = kData{2, 2};

% initialize elementEnergy
elementEnergy.strainEnergy = 0;

%% TANGENT
KWW = KWW + 2 * E * I / he^3 * [6, -6; -6, 6];
KWWApostrophe = KWWApostrophe + 2 * E * I / he^3 * [3 * he, 3 * he; -3 * he, -3 * he];
KWApostropheW = KWApostropheW + 2 * E * I / he^3 * [3 * he, -3 * he; 3 * he, -3 * he];
KWApostropheWApostrophe = KWApostropheWApostrophe + 2 * E * I / he^3 * [2 * he^2, he^2; he^2, 2 * he^2];

%% RESIDUAL
RW = RW + KWW * wN1(:) + KWWApostrophe * wApostropheN1(:);
RWApostrophe = RWApostrophe + KWApostropheW * wN1(:) + KWApostropheWApostrophe * wApostropheN1(:);

%% PASS COMPUTATION DATA
if ~computePostData
    % pass residual
    rData{1} = RW;
    rData{2} = RWApostrophe;

    % pass tangent
    kData{1, 1} = KWW;
    kData{1, 2} = KWWApostrophe;
    kData{2, 1} = KWApostropheW;
    kData{2, 2} = KWApostropheWApostrophe;
end
end
