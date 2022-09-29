function [Ke,Fe] = element(ed,ir,C)
% Elementroutine
%
%   Eingang: ed - Zeile der Elementgeometriedatentabelle
%               ed = [ x1 y1 x2 y2 x3 y3 x4 y4 ];
%            ir - Parameter Integrationsstuetzpunkte
%            C  - Materialmatrix


% Anzahl der Integrationspunkte
ngp = ir*ir;

% Fallunterscheidung - Stuetzpunkte, Gewichte
if (ir == 1)    
    g1 = 0.0; 
    w1 = 2.0;
    gp = [ g1 g1 ];  
    w  = [ w1 w1 ];

elseif (ir == 2)
    g1 = 0.577350269189626; 
    w1 = 1;
    gp(:,1) = [-g1; g1;-g1; g1];  
    gp(:,2) = [-g1;-g1; g1; g1];
    w(:,1)  = [ w1; w1; w1; w1];   
    w(:,2)  = [ w1; w1; w1; w1];

elseif (ir == 3)
    g1 = 0.774596669241483; 
    g2 = 0.;
    w1 = 0.555555555555555; 
    w2 = 0.888888888888888;
    gp(:,1) = [-g1;-g2; g1;-g1; g2; g1;-g1; g2; g1];
    gp(:,2) = [-g1;-g1;-g1; g2; g2; g2; g1; g1; g1];
    w(:,1)  = [ w1; w2; w1; w1; w2; w1; w1; w2; w1];
    w(:,2)  = [ w1; w1; w1; w2; w2; w2; w1; w1; w1];

else
    disp('Used number of integration points not implemented');
    return
end

wp  = w(:,1).*w(:,2);
xi  = gp(:,1);  
eta = gp(:,2);  
r2  = ngp*2;

%% Formfunktionen 
% N = [ N1(xi1,eta1) N2(xi1,eta1) N3(xi1,eta1) N4(xi1,eta1) 
%       N1(xi2,eta2) N2(xi2,eta2) N3(xi2,eta2) N4(xi2,eta2) 
%       N1(xi3,eta3) N2(xi3,eta3) N3(xi3,eta3) N4(xi3,eta3) 
%       N1(xi4,eta4) N2(xi4,eta4) N3(xi4,eta4) N4(xi4,eta4) ]
%
%   Tauchen hier nicht in der STEMA auf, deswegen auskommentiert.

% N(:,1) = (1-xi).*(1-eta)/4;  
% N(:,2) = (1+xi).*(1-eta)/4;
% N(:,3) = (1+xi).*(1+eta)/4;  
% N(:,4) = (1-xi).*(1+eta)/4;     


%% Ableitungen der Formfunktionen
% dNr = [ dN1dxi (eta1)   dN2dxi (eta1)    dN3dxi (eta1)     dN4dxi(eta1)
%         dN1deta(xi1)    dN2deta(xi1)     dN3deta(xi1)      dN4deta(xi1)
%         dN1dxi (eta2)   dN2dxi (eta2)    dN3dxi (eta2)     dN4dxi(eta2)
%         dN1deta(xi2)    dN2deta(xi2)     dN3deta(xi2)      dN4deta(xi2)
%         dN1dxi (eta3)   dN2dxi (eta3)    dN3dxi (eta3)     dN4dxi(eta3)
%         dN1deta(xi3)    dN2deta(xi3)     dN3deta(xi3)      dN4deta(xi3)
%         dN1dxi (eta4)   dN2dxi (eta4)    dN3dxi (eta4)     dN4dxi(eta4)
%         dN1deta(xi4)    dN2deta(xi4)     dN3deta(xi4)      dN4deta(xi4) ]

% Partielle Ableitungen der Formfunktionen nach xi
dNr(1:2:r2,1) = -(1-eta)/4;     
dNr(1:2:r2,2) =  (1-eta)/4;
dNr(1:2:r2,3) =  (1+eta)/4;     
dNr(1:2:r2,4) = -(1+eta)/4;

% Partielle Ableitungen der Formfunktionen nach eta
dNr(2:2:r2+1,1) = -(1-xi)/4;   
dNr(2:2:r2+1,2) = -(1+xi)/4;
dNr(2:2:r2+1,3) =  (1+xi)/4;   
dNr(2:2:r2+1,4) =  (1-xi)/4;

%% Jacobimatrizen (fuer alle Gausspunkte)
%   J = [ J(xi1,eta1) J(xi2,eta2) J(xi3,eta3) J(xi4,eta4) ];
J = [ed(1,(1:2:size(ed,2)-1)); ed(1,(2:2:size(ed,2)))]*dNr';


%% Elementsteifigkeitmatrizen/-lastvektoren aufbauen
% Allozieren
Ke = zeros(8,8);
Fe = zeros(8,1);

for i = 1:ngp
    % Schleife ueber alle Gausspunkte
    
    % Zugriffsindex fuer korrekte Formfunktionen/Jacobimatrix
    indx    = [ 2*i-1; 2*i ];
    
    % Jacobimatrix/-determinante
    JT      = J(:,indx)';
    detJ    = det(JT);

    % Plausibilitaetstest
    if detJ < 10*eps
      disp('Jacobideterminant equal or less than zero!')
    end

    % Ableitung Formfunktion nach globalen Koordinaten
    %   dNdx = [    dN1dx   dN2dx   dN3dx   dN4dx
    %               dN1dy   dN2dy   dN3dy   dN4dy   ]
    dNx = JT\dNr(indx,:);

    % Knotenoperatormatrix
    %   B = [   dN1dx   0       dN2dx   0       dN3dx   0       dN4dx   0
    %           0       dN1dy   0       dN2dy   0       dN3dy   0       dN4dy
    %           dN1dy   dN1dx   dN2dy   dN2dx   dN3dy   dN3dx   dN4dy   dN4dx ]    
    B(1,1:2:8-1) = dNx(1,:);
    B(2,2:2:8)   = dNx(2,:);
    B(3,1:2:8-1) = dNx(2,:);
    B(3,2:2:8)   = dNx(1,:);
    
    % Elementsteifigkeitsmatrix
    Ke = Ke + B'*C*B*detJ*wp(i);
end
