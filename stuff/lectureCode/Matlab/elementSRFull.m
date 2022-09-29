function [KeFull,Fe]=elementSRFull(ed,ir,C)
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
% N(:,1)=(1-xsi).*(1-eta)/4;  
% N(:,2)=(1+xsi).*(1-eta)/4;
% N(:,3)=(1+xsi).*(1+eta)/4;  
% N(:,4)=(1-xsi).*(1+eta)/4;

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

Ke=zeros(8,8);
He=zeros(4,4);
Te=zeros(4,8);
Fe=zeros(12,1);

J011 = J0(1,1); J012 = J0(1,2);
J021 = J0(2,1); J022 = J0(2,2);

F0 = [  J011*J011   J021*J012   2*J011*J012;
        J012*J021   J022*J022   2*J021*J022;
        J011*J021   J012*J022   J011*J022 + J012*J021];
 
j0=det(J0(:,[1 2]));

%% Elementmatrizen
for i=1:ngp
    
    indx = [ 2*i-1; 2*i ];

    JT = J(:,indx)';
    j = det(JT);

    if j<10*eps
      disp('Jacobideterminant equal or less than zero!')
    end
    
    dNx = JT\dNr(indx,:);

    B(1,1:2:8-1)=dNx(1,:);
    B(2,2:2:8)  =dNx(2,:);
    B(3,1:2:8-1)=dNx(2,:);
    B(3,2:2:8)  =dNx(1,:);

    Ke=Ke+B'*C*B*j*wp(i); 
    
    E=[ xsi(i) 0 0 0  ;
        0 eta(i) 0 0 ; 
        0 0 xsi(i) eta(i) ];    
   
    G=j0/j*inv(F0')*E;
    He = He + G'*C*G*j*wp(i);
    Te = Te + G'*C*B*j*wp(i);
end
KeFull = [  Ke  Te';...
            Te  He];