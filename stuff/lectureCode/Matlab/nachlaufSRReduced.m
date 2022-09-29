function [x,y,s,He2] = nachlaufSRReduced(ed,ed_neu,d,ir,C)
%% Anzahl der Integrationspunkte
ngp=ir*ir;

%% Gauss Punkte
if ir==1    
    g1=0.0; 
    w1=2.0;
    gp=[ g1 g1 ];  
    w=[ w1 w1 ];

elseif ir==2
    g1=0.577350269189626; 
    w1=1;
    gp(:,1)=[-g1; g1;-g1; g1];  
    gp(:,2)=[-g1;-g1; g1; g1];
    w(:,1)=[ w1; w1; w1; w1];   
    w(:,2)=[ w1; w1; w1; w1];

elseif ir==3
    g1=0.774596669241483; 
    g2=0.;
    w1=0.555555555555555; 
    w2=0.888888888888888;
    gp(:,1)=[-g1;-g2; g1;-g1; g2; g1;-g1; g2; g1];
    gp(:,2)=[-g1;-g1;-g1; g2; g2; g2; g1; g1; g1];
    w(:,1)=[ w1; w2; w1; w1; w2; w1; w1; w2; w1];
    w(:,2)=[ w1; w1; w1; w2; w2; w2; w1; w1; w1];

else
    disp('Used number of integration points not implemented');
    return
end

wp=w(:,1).*w(:,2);
xsi=gp(:,1);  
eta=gp(:,2);  
r2=ngp*2;

%% Formfunktionen 
N(:,1)=(1-xsi).*(1-eta)/4;  
N(:,2)=(1+xsi).*(1-eta)/4;
N(:,3)=(1+xsi).*(1+eta)/4;  
N(:,4)=(1-xsi).*(1+eta)/4;

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
J = [ed(1,(1:2:8-1));ed(1,(2:2:8))]*dNr';

eta0 = 0;
xsi0 = 0;
% Partielle Ableitungen der Formfunktionen nach xsi
dNr0(1:2:r2,1)=-(1-eta0)/4;     
dNr0(1:2:r2,2)= (1-eta0)/4;
dNr0(1:2:r2,3)= (1+eta0)/4;     
dNr0(1:2:r2,4)=-(1+eta0)/4;
% Partielle Ableitungen der Formfunktionen nach eta
dNr0(2:2:r2+1,1)=-(1-xsi0)/4;   
dNr0(2:2:r2+1,2)=-(1+xsi0)/4;
dNr0(2:2:r2+1,3)= (1+xsi0)/4;   
dNr0(2:2:r2+1,4)= (1-xsi0)/4;
% Jacobimatrix enhanced
J0 = [ed(1,(1:2:8-1));ed(1,(2:2:8))]*dNr0';

J011 = J0(1,1); J012 = J0(1,2);
J021 = J0(2,1); J022 = J0(2,2);

F0 = [  J011*J011   J021*J012   2*J011*J012;
        J012*J021   J022*J022   2*J021*J022;
        J011*J021   J012*J022   J011*J022 + J012*J021];
 
detJ0=det(J0(:,[1 2]));

%% Elementmatrizen
%   Allozieren
He=zeros(4,4);
Te=zeros(4,8);
s  = zeros(4,1);
x  = zeros(4,1);
y  = zeros(4,1);
He2 = zeros(4,4);

% Schleife ueber alle Gausspunkte
for i = 1:ngp

    indx = [ 2*i-1; 2*i ];

    JT = J(:,indx)';
    detJ = det(JT);

    if detJ<10*eps
      disp('Jacobideterminant equal or less than zero!')
    end
    
    dNx = JT\dNr(indx,:);

    B(1,1:2:8-1)=dNx(1,:);
    B(2,2:2:8)  =dNx(2,:);
    B(3,1:2:8-1)=dNx(2,:);
    B(3,2:2:8)  =dNx(1,:);

    E=[ xsi(i)  0       0       0;
        0       eta(i)  0       0;
        0       0       xsi(i)  eta(i)];
    G=detJ0/detJ*inv(F0')*E;
    He = He + G'*C*G*detJ*wp(i);
    Te = Te + G'*C*B*detJ*wp(i);
end

% Schleife ueber alle Gausspunkte
for i = 1:ngp

    indx = [ 2*i-1; 2*i ];

    JT = J(:,indx)';
    detJ = det(JT);

    if detJ<10*eps
      disp('Jacobideterminant equal or less than zero!')
    end
    
    dNx = JT\dNr(indx,:);

    B(1,1:2:8-1)=dNx(1,:);
    B(2,2:2:8)  =dNx(2,:);
    B(3,1:2:8-1)=dNx(2,:);
    B(3,2:2:8)  =dNx(1,:);
    
    E=[ xsi(i)  0       0       0;
        0       eta(i)  0       0;
        0       0       xsi(i)  eta(i)];
    
    G=detJ0/detJ*inv(F0')*E;
        
    % Cauchy Spannungstensor (2d;Voigtsche Notation)
    sigma = C*B*d' + C*G*(-inv(He)*Te*d');
    
    % Vergleichsspannung (Ebener Spannungszustand)
    vonMises = sqrt(sigma(1,1)^2 + sigma(2,1)^2 - sigma(1,1)*sigma(2,1) + 3*sigma(3,1)^2);
    
    s        = s + N(i,:)'*vonMises*detJ*wp(i);
    He2       = He2 + (N(i,:)'*N(i,:))*detJ*wp(i);    
    x(i,1)   = N(i,:)*ed_neu(1,(1:2:size(ed_neu,2)-1))';
    y(i,1)   = N(i,:)*ed_neu(1,(2:2:size(ed_neu,2)))';
end
    