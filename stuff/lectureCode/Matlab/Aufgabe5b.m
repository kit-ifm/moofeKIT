clear all; close all; clc

disp('-------------------------------------------------------------------')
disp('-                        Aufgabe 5b                               -') 
disp('-------------------------------------------------------------------')

%% Konstanten / Vorgaben
%   (Strecken-)Lastvektor
last            = [0;5];                        % last = input('Bitte geben sie den Lastvektor auf den Neumann-Rand an: ');

% Anzahl Elemente
elemente(1,1)   = 10;                           % in x-Richtung
elemente(2,1)   = 10;                           % in y-Richtung

ir              = 2;                            % Anzahl Gausspunkte => ir^2

% Geometrie
hoehe   = 44;                                   % Hoehe
laenge  = 48;                                   % Laenge

% Materialdaten
emod    = 10^2;                                 % E-Modul
nu      = 0.3;               % Querkontraktionszahl (Poissonzahl)
% nu      = 0.4999;

% Materialmatrix
C = emod/((1+nu)*(1-2*nu))*[1-nu   nu   0;
                            nu   1-nu   0;
                            0     0   (1-2*nu)/2];

%% Netzgenerator
% Netz generieren (Elementfreiheitsgradzuordnungstabelle, Geometrievektor)
[edof,q]       = netz(elemente(1,1), elemente(2,1), laenge, hoehe);

% Geometriedatentabellen bzgl. Gebiet (2D Elemente)
[ed]           = extract(edof,q);
% Geometriedatentabellen bzgl. Rand (1D Elemente)
edofRandRechts = edof((end-elemente(2,1)+1):end,3:6);
[edRandRechts] = extract(edofRandRechts,q);

% Ausgangskonfiguration darstellen
draw(ed);

%% Assembly
%   Allozieren / Speicher reservieren
ndof = size(q,1);
K    = sparse(ndof,ndof);
F    = sparse(ndof,1);

% Elementschleife
nel = elemente(1,1)*elemente(2,1);
for i = 1:nel
    [Ke,Fe] = elementSRReduced(ed(i,:),ir,C);
    [K]     = assem(K,Ke,edof(i,:));
    [F]     = assem(F,Fe,edof(i,:));
end

%% Randbedingungen einbauen
% Dirichlet-Rand
RB_1    = 0;
bc      = RB_1*ones(elemente(2,1)*2+2,2);
bc(:,1) = (1:(elemente(2,1)*2+2))';
% Neumann-Rand
RB_2    = last;

% Randelementschleife
for i = 1:elemente(2,1)
    [Fe]    = elementNeumannRand(edRandRechts(i,:),ir,RB_2);
    [F]     = assem(F,Fe,edofRandRechts(i,:));
end

%% Gleichungsloeser
d = solveq(K,F,bc);

%% Darstellung Gleichgewichtslage
q        = q+d;                     % Update, Konfiguration Gleichgewichtslage
[ed_neu] = extract(edof,q);   
draw(ed_neu);

%% Nachlaufrechnung fuer Spannungsberechnung
[d_neu] = extract(edof,d);
s       = sparse(ndof/2,1);
x       = zeros(4,nel);
y       = zeros(4,nel);
edofSpannung = [edof(:,2:2:end)/2];
H            = zeros(ndof/2,ndof/2);
for i = 1:nel
    [x(1:4,i),y(1:4,i),se,He]= nachlaufSRReduced(ed(i,:),ed_neu(i,:),d_neu(i,:),ir,C);
    [s] = assem(s,se,edofSpannung(i,:));
    [H] = assem(H,He,edofSpannung(i,:));
end
spannungVec = H\s;
draw_conf(spannungVec,q,edof,edofSpannung,ed);
colorbar;
colormap jet;