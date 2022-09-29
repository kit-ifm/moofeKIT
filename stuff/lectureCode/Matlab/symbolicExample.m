syms X1 X2 X3

%displacement (analytic)
xSymbolic(X1,X2,X3) = [X1+10*X1^2; X2; X3];
    
Fana = GradientOfVec(xSymbolic);
Fana = Fana(X1,X2,X3);
% constitutive equations
% Psi = a*(trace(F'*F)-3);
DPsi_F = 12*2*Fana;
Pana = DPsi_F;
% body force
BSymbolic = -DivergenceOfTensor(Pana);

%% SUBFUNCTIONS
function out = GradientOfVec(u)
    syms X1 X2 X3 
    X = [X1;X2;X3];
    out = symfun(zeros(3,3),[X1,X2,X3]);
    for i = 1:3
        out = out + horzcat(zeros(3,i-1),diff(u,X(i)),zeros(3,3-i));
    end
end

function out = DivergenceOfTensor(A)
    syms X1 X2 X3 
    X = [X1;X2;X3];
    out = symfun(zeros(3,1),[X1,X2,X3]);
    I = eye(3);
    for i=1:3
        out = out + diff(A,X(i))* I(:,i);
    end
end
