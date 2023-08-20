function [x,y,s,He] = nachlauf(ed,ed_neu,d,ir,C)


% Anzahl der Integrationspunkte
ngp=ir*ir;

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

wp = w(:,1).*w(:,2);

xsi = gp(:,1);  
eta = gp(:,2);  
r2  = ngp*2;

%% Formfunktionen 
% N = [ N1(xi1,eta1) N2(xi1,eta1) N3(xi1,eta1) N4(xi1,eta1) 
%       N1(xi2,eta2) N2(xi2,eta2) N3(xi2,eta2) N4(xi2,eta2) 
%       N1(xi3,eta3) N2(xi3,eta3) N3(xi3,eta3) N4(xi3,eta3) 
%       N1(xi4,eta4) N2(xi4,eta4) N3(xi4,eta4) N4(xi4,eta4) ]
%
N(:,1)=(1-xsi).*(1-eta)/4;  
N(:,2)=(1+xsi).*(1-eta)/4;
N(:,3)=(1+xsi).*(1+eta)/4;  
N(:,4)=(1-xsi).*(1+eta)/4;

%% Ableitungen der Formfunktionen
% dNr = [ dN1dxi (eta1)   dN2dxi (eta1)    dN3dxi (eta1)     dN4dxi(eta1)
%         dN1deta(xi1)    dN2deta(xi1)     dN3deta(xi1)      dN4deta(xi1)
%         dN1dxi (eta2)   dN2dxi (eta2)    dN3dxi (eta2)     dN4dxi(eta2)
%         dN1deta(xi2)    dN2deta(xi2)     dN3deta(xi2)      dN4deta(xi2)
%         dN1dxi (eta3)   dN2dxi (eta3)    dN3dxi (eta3)     dN4dxi(eta3)
%         dN1deta(xi3)    dN2deta(xi3)     dN3deta(xi3)      dN4deta(xi3)
%         dN1dxi (eta4)   dN2dxi (eta4)    dN3dxi (eta4)     dN4dxi(eta4)
%         dN1deta(xi4)    dN2deta(xi4)     dN3deta(xi4)      dN4deta(xi4) ]

% Partielle Ableitungen der Formfunktionen nach xsi
dNr(1:2:r2,1)=-(1-eta)/4;     
dNr(1:2:r2,2)= (1-eta)/4;
dNr(1:2:r2,3)= (1+eta)/4;     
dNr(1:2:r2,4)=-(1+eta)/4;

% Partielle Ableitungen der Formfunktionen nach eta
dNr(2:2:r2+1,1)=-(1-xsi)/4;   
dNr(2:2:r2+1,2)=-(1+xsi)/4;
dNr(2:2:r2+1,3)= (1+xsi)/4;   
dNr(2:2:r2+1,4)= (1-xsi)/4;

% Jacobimatrix
%   J = [ J(xi1,eta1) J(xi2,eta2) J(xi3,eta3) J(xi4,eta4) ];
% Ausgangskonfiguration
J = [ed(1,(1:2:size(ed,2)-1));ed(1,(2:2:size(ed,2)))]*dNr';

%% Elementmatrizen
%   Allozieren
s  = zeros(4,1);
x  = zeros(4,1);
y  = zeros(4,1);
He = zeros(4,4);

for i = 1:ngp
    % Schleife ueber alle Gausspunkte

    % Zugriffsindex fuer korrekte Formfunktionen/Jacobimatrix
    indx = [ 2*i-1; 2*i ];
    
    % Jakobimatrizen/-determinante
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
    
    % Cauchy Spannungstensor (2d;Voigtsche Notation)
    sigma   = C*B*d';
    
    % Vergleichsspannung (Ebener Spannungszustand)
    vonMises = sqrt(sigma(1,1)^2 + sigma(2,1)^2 - sigma(1,1)*sigma(2,1) + 3*sigma(3,1)^2);
    
    s        = s + N(i,:)'*vonMises*detJ*wp(i);
    He       = He + (N(i,:)'*N(i,:))*detJ*wp(i);    
    x(i,1)   = N(i,:)*ed_neu(1,(1:2:size(ed_neu,2)-1))';
    y(i,1)   = N(i,:)*ed_neu(1,(2:2:size(ed_neu,2)))';
end
