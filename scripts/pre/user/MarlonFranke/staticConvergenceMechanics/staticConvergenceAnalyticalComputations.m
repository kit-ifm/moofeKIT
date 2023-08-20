function [x,C,G,c,lambdaC,lambdaG,lambdac,B] = staticConvergenceAnalyticalComputations(a,b,c1,d1,g1,g2,g3,choiceFunction,dimension)
%% calculation of analytic bodyforce
syms X1 X2 X3

%displacement (analytic)
if choiceFunction == 1
    x(X1,X2,X3) = [X1+g1*X1^2; X2; X3];                              
elseif choiceFunction == 2
    x(X1,X2,X3) = [X1+g1*sin(X1); X2+g2*cos(X2); X3+g3*(sin(X3)+cos(X3))];   %from Poya et al.
elseif choiceFunction == 3
    x(X1,X2,X3) = [X1+g1*X1^3;X2+g2*X2^3;X3+g3*X3^3];                        %from Ortigosa et al. 2016 
elseif choiceFunction == 4
    x(X1,X2,X3) = [X1+g1*X1^2; X2+g2*X2^2; X3+g3*X3^2];
elseif choiceFunction == 5
    x(X1,X2,X3) = [X1;X2;X3];
end

% kinematics
F = GradientOfVector(x);
F = F(X1,X2,X3);
% Hana = 1/2*wedge(Fana,Fana);
% Jana = det(Fana);

C = F.'*F;
G = 1/2*wedge(C,C);
c = 1/3*sum(sum(C.*G));

% constitutive equations
%     Psi = a*(trace(Cana)-3) + b*(trace(Gana)-3) + c1/2*(sqrt(cana)-1)^2 - d1*log(sqrt(cana));
%     Sana = 2*(a*eye(3) + 1/(2*er*e0*cana^(1/2))*(Dana*Dana') + wedge(b*eye(3),Cana) + (c1*(1-1/sqrt(cana)) - d1/(2*cana) - 1/(4*er*e0)*cana^(-3/2)*Dana'*(Cana*Dana))*Gana);
% strainEnergy = (a * (trace(CxN1) - 3) + b * (trace(GxN1) - 3) - d1 * log(sqrt(cxN1)) + c1 / 2 * (sqrt(cxN1) - 1)^2);
DW_C = a*eye(3);
DW_G = b*eye(3);
DW_c = c1/2*(1-1/sqrt(c)) - d1/(2*c);

% TODO: why formula?
lambdac = formula(DW_c);
lambdaG = DW_G + 1/3*lambdac*C;
% TODO 1. formula, 2. DW_C
lambdaC = DW_C + wedge(lambdaG,C) + 1/3*lambdac*G;

S = 2*lambdaC;

P(X1,X2,X3) = F*S;
% body force
B = -DivergenceOfTensor(P);

P2(X1,X2,X3) = F*S;
xtest = rand(3,1);
PNumerical = matlabFunction(P);
P2Numerical = matlabFunction(P2);
toleranceAnalyticalStressCalculation = 9;
if any(any(round(PNumerical(xtest(1),xtest(2),xtest(3))-P2Numerical(xtest(1),xtest(2),xtest(3)),toleranceAnalyticalStressCalculation)~=0)~=0)
    disp(PNumerical(xtest(1),xtest(2),xtest(3))-P2Numerical(xtest(1),xtest(2),xtest(3)))
    error('stress tensors wrong')
else
    clearvars xtest PanaNum Pana2Num
end
% symbolic2Function
x = matlabFunction(x);
C(X1,X2,X3) = C;
C = matlabFunction(C);
G(X1,X2,X3) = G;
G = matlabFunction(G);
c(X1,X2,X3) = c;
c = matlabFunction(c);
lambdaC(X1,X2,X3) = lambdaC;
lambdaC = matlabFunction(lambdaC);
lambdaG(X1,X2,X3) = lambdaG;
lambdaG = matlabFunction(lambdaG);
lambdac(X1,X2,X3) = lambdac;
lambdac = matlabFunction(lambdac);
B = matlabFunction(B);
end

%% SUBFUNCTIONS
function out=GradientOfVector(u)
    syms X1 X2 X3 
    X = [X1;X2;X3];
    out=symfun(zeros(3,3),[X1,X2,X3]);
    for i=1:3
        out = out + horzcat(zeros(3,i-1),diff(u,X(i)),zeros(3,3-i));
    end
end

function out=GradientOfScalar(u)
    syms X1 X2 X3 
    X = [X1;X2;X3];
    out=symfun(zeros(3,1),[X1,X2,X3]);
    I = eye(3);
    for i=1:3
        out = out + diff(u,X(i))* I(:,i);
    end
end

function out=DivergenceOfTensor(A)
    syms X1 X2 X3 
    X = [X1;X2;X3];
    out=symfun(zeros(3,1),[X1,X2,X3]);
    I = eye(3);

    for i=1:3
        out = out + diff(A,X(i))* I(:,i);
    end
end

function out=DivergenceOfVector(A)
    syms X1 X2 X3 
    X = [X1;X2;X3];
    out=symfun(zeros(1,1),[X1,X2,X3]);
    I = eye(3);

    for i=1:3
        out = out + diff(A,X(i))'*I(:,i);
    end
end

function [H] = wedge(A,B)  
%% Creates tensor cross product.
    %
    % Syntax
    % A,B vector or matrix
    % [H] = wedge(A,B)      
    %
    % 22.11.2017 Robin Pfefferkorn
    
    H = zeros(size(A));
    symbolic = ~isempty(symvar(A(:)')) || ~isempty(symvar(B(:)'));
    if symbolic
        H = sym(H);
    end
    
    % assertions of size
    assert(size(A,1)==3 && (size(A,2)==1 || size(A,2)==3),'A must have 3 rows and either 1 or 3 columns')
    assert(size(B,1)==3 && (size(B,2)==1 || size(B,2)==3),'B must have 3 rows and either 1 or 3 columns')
    assert(size(A,2)==3 || size(B,2)==3,'A or B must be a 3x3 matrix')
   
    if size(A,2)==3 && size(B,2)==3
        H(1,:) = cross(A(2,:),B(3,:)) + cross(B(2,:),A(3,:));
        H(2,:) = cross(B(3,:),A(1,:)) + cross(A(3,:),B(1,:));
        H(3,:) = cross(A(1,:),B(2,:)) + cross(B(1,:),A(2,:));
    elseif size(A,2)==1 && size(B,2)==3
        H(:,1) = cross(A,B(:,1));
        H(:,2) = cross(A,B(:,2));
        H(:,3) = cross(A,B(:,3));
    elseif size(A,2)==3 && size(B,2)==1
        H(1,:) = cross(A(1,:),B);
        H(2,:) = cross(A(2,:),B);
        H(3,:) = cross(A(3,:),B);
    else
        error('Matrix dimensions do not agree!')
    end
     
end
