function [Fe]=element_rand(ed_rand,ir,RB_2)
% Lastroutine (Neumannrand)
%
%   siehe Elementroutine


% Anzahl der Integrationspunkte
ngp = ir;

% Fallunterscheidung - Stuetzpunkte, Gewichte
if ir == 1
    g1 = 0.0; 
    w1 = 2.0;
    gp = g1; 
    w  = w1;
    
elseif ir == 2
    g1 = 0.577350269189626; 
    w1 = 1;
    gp = [-g1; g1]; 
    w  = [ w1; w1];    
        
elseif ir == 3
    g1 = 0.774596669241483; 
    g2 = 0;
    w1 = 0.555555555555555; 
    w2 = 0.888888888888888;
    gp = [-g1; g2; g1];
    w  = [ w1; w2; w1];
   
else
    disp('Used number of intergratins poits not implemented')
    return
end

xsi = gp(:);  

%% Lineare Formfunktionen
%   N = [   N1(xi1)     N2(xi1)
%           N1(xi2)     N2(xi2) ]
N(:,1)=(1-xsi)/2;  
N(:,2)=(1+xsi)/2;
    
% Elementlaenge
Le = norm([ed_rand(3);ed_rand(4)] - [ed_rand(1);ed_rand(2)]);

% Jacobimatrix/-determinante
J = Le/2;
detJ = det(J);
    
%% Elementlastvektor aufbauen
Fe = zeros(4,1);
for ii = 1:ngp
    % Schleife ueber alle Gausspunkte   
    Fe([1;3]) = Fe([1;3]) + N(ii,:)'*RB_2(1)*detJ*w(ii);
    Fe([2;4]) = Fe([2;4]) + N(ii,:)'*RB_2(2)*detJ*w(ii);
end
    
   
        
