function [O] = assem(O,Oe,dof)
%% Finite element assembly
% 
% K = assem(K,Ke,dof)
% F = assem(F,Fe,dof)
%
% K: Steifigkeitsmatrix
% Ke: Elementsteifigkeitsmatrix
%
% F: Lastvektor
% Fe: Elementlastvektor
%
% dof: Indices der Elementfreiheitsgrade
%
% Zusammenbau der Finite-Elemente-Matrizes oder FE-Vektoren

[nRows,nCols] = size(Oe);

if nRows == nCols
    % Matrix-Assembly
    % K(dof,dof) = K(dof,dof) + Ke
    O(dof,dof) = O(dof,dof) + Oe;
else
    % Vektor-Assembly
    % F(dof) = F(dof) + Fe;
    O(dof) = O(dof) + Oe;
end

end