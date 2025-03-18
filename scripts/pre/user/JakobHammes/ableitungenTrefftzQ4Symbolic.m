% Ableitungen der pQ-Matrix mittels MATLAB
syms x y aXi aEta bXi bEta x0 y0 nuBar 

%Definition x/y-Dach
xHat = x-x0;
yHat = y-y0;

% Definition "Hilfskoordinatensystem" 
rXi  = (1/sqrt(aXi^2 + bXi^2))*(aXi*xHat+bXi*yHat);
sXi  = (1/sqrt(aXi^2 + bXi^2))*(-bXi*xHat+aXi*yHat);

rEta = (1/sqrt(aEta^2 + bEta^2))*(aEta*xHat+bEta*yHat);
sEta = (1/sqrt(aEta^2 + bEta^2))*(-bEta*xHat+aEta*yHat);

% Konstante Matrizen
tXi  = 1/sqrt(aXi^2 + bXi^2)*[aXi -bXi; bXi aXi];
tEta = 1/sqrt(aEta^2 + bEta^2)*[aEta -bEta; bEta aEta];

% pQ-Matrix
pQ = [1 0 xHat 0 yHat 0 tXi(1,:)*[-2*rXi*sXi; nuBar*sXi^2+rXi^2] tEta(1,:)*[-2*rEta*sEta; nuBar*sEta^2+rEta^2]
      0 1 0 xHat 0 yHat tXi(2,:)*[-2*rXi*sXi; nuBar*sXi^2+rXi^2] tEta(2,:)*[-2*rEta*sEta; nuBar*sEta^2+rEta^2]];


%Ableiten nach x/y
pQdx = diff(pQ,x);
pQdy = diff(pQ,y);


%Ausgabe
disp('dpQ/dx')
disp(pQdx)
disp('-------------------------------')
disp('dpQ/dy')
disp(pQdy)


% Ergebnisse für Elementformulierung aufbereitet 
% [0, 0, 1, 0, 0, 0, (2*aXi^2*(bXi*(xHat) - aXi*(yHat)))/(aXi^2 + bXi^2)^(3/2) - (bXi*((2*aXi*(aXi*(xHat) + bXi*(yHat)))/(aXi^2 + bXi^2) + (2*bXi*nuBar*(bXi*(xHat) - aXi*(yHat)))/(aXi^2 + bXi^2)))/(aXi^2 + bXi^2)^(1/2) + (2*aXi*bXi*(aXi*(xHat) + bXi*(yHat)))/(aXi^2 + bXi^2)^(3/2), (2*aEta^2*(bEta*(xHat) - aEta*(yHat)))/(aEta^2 + bEta^2)^(3/2) - (bEta*((2*aEta*(aEta*(xHat) + bEta*(yHat)))/(aEta^2 + bEta^2) + (2*bEta*nuBar*(bEta*(xHat) - aEta*(yHat)))/(aEta^2 + bEta^2)))/(aEta^2 + bEta^2)^(1/2) + (2*aEta*bEta*(aEta*(xHat) + bEta*(yHat)))/(aEta^2 + bEta^2)^(3/2)
%  0, 0, 0, 1, 0, 0, (aXi*((2*aXi*(aXi*(xHat) + bXi*(yHat)))/(aXi^2 + bXi^2) + (2*bXi*nuBar*(bXi*(xHat) - aXi*(yHat)))/(aXi^2 + bXi^2)))/(aXi^2 + bXi^2)^(1/2) + (2*bXi^2*(aXi*(xHat) + bXi*(yHat)))/(aXi^2 + bXi^2)^(3/2) + (2*aXi*bXi*(bXi*(xHat) - aXi*(yHat)))/(aXi^2 + bXi^2)^(3/2), (aEta*((2*aEta*(aEta*(xHat) + bEta*(yHat)))/(aEta^2 + bEta^2) + (2*bEta*nuBar*(bEta*(xHat) - aEta*(yHat)))/(aEta^2 + bEta^2)))/(aEta^2 + bEta^2)^(1/2) + (2*bEta^2*(aEta*(xHat) + bEta*(yHat)))/(aEta^2 + bEta^2)^(3/2) + (2*aEta*bEta*(bEta*(xHat) - aEta*(yHat)))/(aEta^2 + bEta^2)^(3/2)]
 

%[0, 0, 0, 0, 1, 0, (2*aXi*bXi*(bXi*(xHat) - aXi*(yHat)))/(aXi^2 + bXi^2)^(3/2) - (2*aXi^2*(aXi*(xHat) + bXi*(yHat)))/(aXi^2 + bXi^2)^(3/2) - (bXi*((2*bXi*(aXi*(xHat) + bXi*(yHat)))/(aXi^2 + bXi^2) - (2*aXi*nuBar*(bXi*(xHat) - aXi*(yHat)))/(aXi^2 + bXi^2)))/(aXi^2 + bXi^2)^(1/2), (2*aEta*bEta*(bEta*(xHat) - aEta*(yHat)))/(aEta^2 + bEta^2)^(3/2) - (2*aEta^2*(aEta*(xHat) + bEta*(yHat)))/(aEta^2 + bEta^2)^(3/2) - (bEta*((2*bEta*(aEta*(xHat) + bEta*(yHat)))/(aEta^2 + bEta^2) - (2*aEta*nuBar*(bEta*(xHat) - aEta*(yHat)))/(aEta^2 + bEta^2)))/(aEta^2 + bEta^2)^(1/2)
% 0, 0, 0, 0, 0, 1, (aXi*((2*bXi*(aXi*(xHat) + bXi*(yHat)))/(aXi^2 + bXi^2) - (2*aXi*nuBar*(bXi*(xHat) - aXi*(yHat)))/(aXi^2 + bXi^2)))/(aXi^2 + bXi^2)^(1/2) + (2*bXi^2*(bXi*(xHat) - aXi*(yHat)))/(aXi^2 + bXi^2)^(3/2) - (2*aXi*bXi*(aXi*(xHat) + bXi*(yHat)))/(aXi^2 + bXi^2)^(3/2), (aEta*((2*bEta*(aEta*(xHat) + bEta*(yHat)))/(aEta^2 + bEta^2) - (2*aEta*nuBar*(bEta*(xHat) - aEta*(yHat)))/(aEta^2 + bEta^2)))/(aEta^2 + bEta^2)^(1/2) + (2*bEta^2*(bEta*(xHat) - aEta*(yHat)))/(aEta^2 + bEta^2)^(3/2) - (2*aEta*bEta*(aEta*(xHat) + bEta*(yHat)))/(aEta^2 + bEta^2)^(3/2)]


