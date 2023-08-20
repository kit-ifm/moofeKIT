function q = solveq(K,F,bc)
%% Loesen des lineare Gleichungssystems unter Beruecksichtiung Dirichlet-Randbedingungen
%
% q = solveq(K,F,bc)
%
% q: Loesungsvektor
% K: Steifigkeitsmatrix
% F: Lastvektor
% bc: Dirichlet-Freiheitsgrade
%
% bc = [idxDiri1 valueDiri1
%       idxDiri2 valueDiri2
%            :              
%       idxDiriN valueDiriN]

if ~isempty(bc)
    % Eingaben einlesen
    dofDiri = bc(:,1);
    qDiri   = bc(:,2);

    % Freiheitsgrade mit Dirichlet-Raender identifizieren
    nDof = numel(F);
    DOSOLVE = true(nDof,1);
    DOSOLVE(dofDiri) = false;
    
    % Nach Freiheitsgraden, die nicht Dirichlet-Rand sind loesen, bei
    % anderen Werten Dirichletrand einsetzen
    q = zeros(size(K,1),1);
    q(DOSOLVE) = K(DOSOLVE,DOSOLVE)\(F(DOSOLVE) - K(DOSOLVE,~DOSOLVE)*qDiri);
    q(~DOSOLVE) = qDiri;
else
    % Nach allen Werten loesen
    q = K\F;
end

end