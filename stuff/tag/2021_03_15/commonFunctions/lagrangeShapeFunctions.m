function out = lagrangeShapeFunctions(numNodes,numberOfGausspoints,dimension)
%ansatzfunctions at gausspoints
out = struct('N',[],'dNr',[],'dNr2',[],'gaussWeight',[],'gaussPoint',[]);

g1 = 0.577350269189626;
g2 = 0.774596669241483;
g3 = 0.339981043584856;
g4 = 0.861136311594053;
     
w1 = 0.555555555555556;
w2 = 0.888888888888889;
w3 = 0.652145154862546;
w4 = 0.347854845137454;

if dimension == 1
    % Gausspoints
    if numberOfGausspoints == 1
        gaussPoint = 0;
        gaussWeight = 2;
    elseif numberOfGausspoints == 2
        g = 1/sqrt(3);
        gaussPoint = [-g, g];
        gaussWeight = [1, 1];
    elseif numberOfGausspoints == 3
        g = sqrt(3/5);
        gaussPoint = [-g, 0, g];
        gaussWeight = [5/9,8/9,5/9];
    else
        error('number of Gausspoints not implemented')
    end

    % Ansatzfunctions at gausspoints
    xsi = gaussPoint(1,:);
    N = zeros(numberOfGausspoints,numNodes);
    dNr = N;
    dNr2 = N;
    
    if numNodes == 1
        N(:,1) = 1;
        dNr(:,1) = 0;
        dNr2(:,1) = 0;
    elseif numNodes == 2
        % shape functions
        N(:,1)=1/2*(1-xsi);
        N(:,2)=1/2*(1+xsi);
        % derivativ wrt xsi
        dNr(:,1)=-1/2;
        dNr(:,2)= 1/2;
        % second derivativ wrt xsi
        dNr2(:,1) = 0;
        dNr2(:,2) = 0;
    elseif numNodes == 3
        % shape functions
        N(:,1)=1/2*xsi.*(xsi-1);
        N(:,2)=1/2*xsi.*(xsi+1);
        N(:,3)=1-xsi.^2;
        % derivativ wrt xsi
        dNr(:,1) = xsi-1/2;
        dNr(:,2) = xsi+1/2;
        dNr(:,3) = -2*xsi;
        % second derivativ wrt xsi
        dNr2(:,1) = 1;
        dNr2(:,2) = 1;
        dNr2(:,3) = -2;
    else
        error('Order of ansatzfunctions not implemented'); 
    end
elseif dimension == 2
    % Gausspoints
    if numberOfGausspoints == 1
        gaussPoint = [0; 0];
        gaussWeight = 2;
    elseif numberOfGausspoints == 4
        g = 1/sqrt(3);
        gaussPoint = [-g, g, g,-g;
            -g,-g, g, g];
        gaussWeight = [1, 1, 1, 1];
    elseif numberOfGausspoints == 9
        g = sqrt(3/5);
        gaussPoint = [-g, 0, g,-g, 0, g,-g, 0, g;
            -g,-g,-g, 0, 0, 0, g, g, g];
        gaussWeight = [25/81,40/81,25/81,40/81,64/81,40/81,25/81,40/81,25/81];
    else
        error('number of Gausspoints not implemented')
    end
    
    % Ansatzfunctions at gausspoints
    xsi = gaussPoint(1,:);
    eta = gaussPoint(2,:);
    r2 = numberOfGausspoints*2;
    r3 = numberOfGausspoints*3;
    N = zeros(numberOfGausspoints,numNodes);
    dNr = zeros(r2,numNodes);
    dNr2 = zeros(r3,numNodes);
    
    if numNodes == 1
        N(:,:) = 1;
        dNr(:,:) = 0;
        dNr2(:,:) = 0;
    elseif numNodes == 4
        % shape functions
        N(:,1)=(1-xsi).*(1-eta)/4;
        N(:,2)=(1+xsi).*(1-eta)/4;
        N(:,3)=(1+xsi).*(1+eta)/4;
        N(:,4)=(1-xsi).*(1+eta)/4;
        % derivativ wrt xsi
        dNr(1:2:r2,1)=-(1-eta)/4;
        dNr(1:2:r2,2)= (1-eta)/4;
        dNr(1:2:r2,3)= (1+eta)/4;
        dNr(1:2:r2,4)=-(1+eta)/4;
        % derivative wrt eta
        dNr(2:2:r2,1)=-(1-xsi)/4;
        dNr(2:2:r2,2)=-(1+xsi)/4;
        dNr(2:2:r2,3)= (1+xsi)/4;
        dNr(2:2:r2,4)= (1-xsi)/4;
    elseif numNodes == 9
        % shape functions
        N(:,1)=1/4*(xsi.^2-xsi).*(eta.^2-eta);
        N(:,2)=1/4*(xsi.^2+xsi).*(eta.^2-eta);
        N(:,3)=1/4*(xsi.^2+xsi).*(eta.^2+eta);
        N(:,4)=1/4*(xsi.^2-xsi).*(eta.^2+eta);
        N(:,5)=1/2*(1-xsi.^2).*eta.*(eta-1);
        N(:,6)=1/2*xsi.*(xsi+1).*(1-eta.^2);
        N(:,7)=1/2*(1-xsi.^2).*eta.*(eta+1);
        N(:,8)=1/2*xsi.*(xsi-1).*(1-eta.^2);
        N(:,9)=(1-xsi.^2).*(1-eta.^2);
        % derivativ wrt xsi
        dNr(1:2:r2,1) = 1/4*(2*xsi-1).*(eta.^2-eta);
        dNr(1:2:r2,2) = 1/4*(2*xsi+1).*(eta.^2-eta);
        dNr(1:2:r2,3) = 1/4*(2*xsi+1).*(eta.^2+eta);
        dNr(1:2:r2,4) = 1/4*(2*xsi-1).*(eta.^2+eta);
        dNr(1:2:r2,5) = -xsi.*eta.*(eta-1);
        dNr(1:2:r2,6) = 1/2*(2*xsi+1).*(1-eta.^2);
        dNr(1:2:r2,7) = -xsi.*eta.*(eta+1);
        dNr(1:2:r2,8) = 1/2*(2*xsi-1).*(1-eta.^2);
        dNr(1:2:r2,9) = -2*xsi.*(1-eta.^2);
        % derivative wrt eta
        dNr(2:2:r2,1) = 1/4*(xsi.^2-xsi).*(2*eta-1);
        dNr(2:2:r2,2) = 1/4*(xsi.^2+xsi).*(2*eta-1);
        dNr(2:2:r2,3) = 1/4*(xsi.^2+xsi).*(2*eta+1);
        dNr(2:2:r2,4) = 1/4*(xsi.^2-xsi).*(2*eta+1);
        dNr(2:2:r2,5) = 1/2*(1-xsi.^2).*(2*eta-1);
        dNr(2:2:r2,6) = -xsi.*(xsi+1).*eta;
        dNr(2:2:r2,7) = 1/2*(1-xsi.^2).*(2*eta+1);
        dNr(2:2:r2,8) = -xsi.*(xsi-1).*eta;
        dNr(2:2:r2,9) = -2*(1-xsi.^2).*eta;
    else
        error('Order of ansatzfunctions not implemented'); 
    end
elseif dimension == 3
    if numberOfGausspoints == 1
        gaussPoint = [0 0 0];
        gaussWeight = 8;
    elseif numberOfGausspoints == 8
        gaussPoint=[-g1 -g1 -g1;
            g1 -g1 -g1;
            -g1  g1 -g1;
            g1  g1 -g1;
            -g1 -g1  g1;
            g1 -g1  g1;
            -g1  g1  g1;
            g1  g1  g1];
        gaussWeight = ones(8,1);
    elseif numberOfGausspoints == 27
        gaussPoint =[-g2 0 g2 -g2 0 g2 -g2 0 g2 -g2 0 g2 -g2 0 g2 -g2 0 g2 -g2 0 g2 -g2 0 g2 -g2 0 g2;
            -g2 -g2 -g2 0 0 0 g2 g2 g2 -g2 -g2 -g2 0 0 0 g2 g2 g2 -g2 -g2 -g2 0 0 0 g2 g2 g2;
            -g2 -g2 -g2 -g2 -g2 -g2 -g2 -g2 -g2 0 0 0 0 0 0 0 0 0 g2 g2 g2 g2 g2 g2 g2 g2 g2]';
        gaussWeight = [[25/81 40/81 25/81 40/81 64/81 40/81 25/81 40/81 25/81]'*w1;
            [25/81 40/81 25/81 40/81 64/81 40/81 25/81 40/81 25/81]'*w2;
            [25/81 40/81 25/81 40/81 64/81 40/81 25/81 40/81 25/81]'*w1];
    else
        error('Number of Gauss points not implemented');
    end
    % Support
    r3 = numberOfGausspoints*3;
    xsi=gaussPoint(:,1);
    eta=gaussPoint(:,2);
    zeta=gaussPoint(:,3);
    if numNodes == 8
        % Shape functions
        N(:,1)=(1-xsi).*(1-eta).*(1-zeta)/8;
        N(:,2)=(1+xsi).*(1-eta).*(1-zeta)/8;
        N(:,3)=(1+xsi).*(1+eta).*(1-zeta)/8;
        N(:,4)=(1-xsi).*(1+eta).*(1-zeta)/8;
        N(:,5)=(1-xsi).*(1-eta).*(1+zeta)/8;
        N(:,6)=(1+xsi).*(1-eta).*(1+zeta)/8;
        N(:,7)=(1+xsi).*(1+eta).*(1+zeta)/8;
        N(:,8)=(1-xsi).*(1+eta).*(1+zeta)/8;
        
        % Derivativ wrt xsi
        dNr(1:3:r3,1)=-(1-eta).*(1-zeta)/8;
        dNr(1:3:r3,2)= (1-eta).*(1-zeta)/8;
        dNr(1:3:r3,3)= (1+eta).*(1-zeta)/8;
        dNr(1:3:r3,4)=-(1+eta).*(1-zeta)/8;
        dNr(1:3:r3,5)=-(1-eta).*(1+zeta)/8;
        dNr(1:3:r3,6)= (1-eta).*(1+zeta)/8;
        dNr(1:3:r3,7)= (1+eta).*(1+zeta)/8;
        dNr(1:3:r3,8)=-(1+eta).*(1+zeta)/8;
        
        % Derivative wrt eta
        dNr(2:3:r3,1)=-(1-xsi).*(1-zeta)/8;
        dNr(2:3:r3,2)=-(1+xsi).*(1-zeta)/8;
        dNr(2:3:r3,3)= (1+xsi).*(1-zeta)/8;
        dNr(2:3:r3,4)= (1-xsi).*(1-zeta)/8;
        dNr(2:3:r3,5)=-(1-xsi).*(1+zeta)/8;
        dNr(2:3:r3,6)=-(1+xsi).*(1+zeta)/8;
        dNr(2:3:r3,7)= (1+xsi).*(1+zeta)/8;
        dNr(2:3:r3,8)= (1-xsi).*(1+zeta)/8;
        
        % Derivative wrt zeta
        dNr(3:3:r3,1)=-(1-xsi).*(1-eta)/8;
        dNr(3:3:r3,2)=-(1+xsi).*(1-eta)/8;
        dNr(3:3:r3,3)=-(1+xsi).*(1+eta)/8;
        dNr(3:3:r3,4)=-(1-xsi).*(1+eta)/8;
        dNr(3:3:r3,5)= (1-xsi).*(1-eta)/8;
        dNr(3:3:r3,6)= (1+xsi).*(1-eta)/8;
        dNr(3:3:r3,7)= (1+xsi).*(1+eta)/8;
        dNr(3:3:r3,8)= (1-xsi).*(1+eta)/8;
        
        dNr2(1,1:3:r3,1) = 0;
        dNr2(1,1:3:r3,2) = 0;
        dNr2(1,1:3:r3,3) = 0;
        dNr2(1,1:3:r3,4) = 0;
        dNr2(1,1:3:r3,5) = 0;
        dNr2(1,1:3:r3,6) = 0;
        dNr2(1,1:3:r3,7) = 0;
        dNr2(1,1:3:r3,8) = 0;
        
        dNr2(1,2:3:r3,1) = (1-zeta)/8;
        dNr2(1,2:3:r3,2) = -(1-zeta)/8;
        dNr2(1,2:3:r3,3) = (1-zeta)/8;
        dNr2(1,2:3:r3,4) = -(1-zeta)/8;
        dNr2(1,2:3:r3,5) = (1+zeta)/8;
        dNr2(1,2:3:r3,6) = -(1+zeta)/8;
        dNr2(1,2:3:r3,7) = (1+zeta)/8;
        dNr2(1,2:3:r3,8) = -(1+zeta)/8;
        
        dNr2(1,3:3:r3,1) = (1-eta)/8;
        dNr2(1,3:3:r3,2) = -(1-eta)/8;
        dNr2(1,3:3:r3,3) = -(1+eta)/8;
        dNr2(1,3:3:r3,4) = (1+eta)/8;
        dNr2(1,3:3:r3,5) = -(1-eta)/8;
        dNr2(1,3:3:r3,6) = (1-eta)/8;
        dNr2(1,3:3:r3,7) = (1+eta)/8;
        dNr2(1,3:3:r3,8) = -(1+eta)/8;
        
        dNr2(2,1:3:r3,1) = (1-zeta)/8;
        dNr2(2,1:3:r3,2) = -(1-zeta)/8;
        dNr2(2,1:3:r3,3) = (1-zeta)/8;
        dNr2(2,1:3:r3,4) = -(1-zeta)/8;
        dNr2(2,1:3:r3,5) = (1+zeta)/8;
        dNr2(2,1:3:r3,6) = -(1+zeta)/8;
        dNr2(2,1:3:r3,7) = (1+zeta)/8;
        dNr2(2,1:3:r3,8) = -(1+zeta)/8;
        
        dNr2(2,2:3:r3,1) = 0;
        dNr2(2,2:3:r3,2) = 0;
        dNr2(2,2:3:r3,3) = 0;
        dNr2(2,2:3:r3,4) = 0;
        dNr2(2,2:3:r3,5) = 0;
        dNr2(2,2:3:r3,6) = 0;
        dNr2(2,2:3:r3,7) = 0;
        dNr2(2,2:3:r3,8) = 0;
        
        dNr2(2,3:3:r3,1) = (1-xsi)/8;
        dNr2(2,3:3:r3,2) = (1+xsi)/8;
        dNr2(2,3:3:r3,3) = -(1+xsi)/8;
        dNr2(2,3:3:r3,4) = -(1-xsi)/8;
        dNr2(2,3:3:r3,5) = -(1-xsi)/8;
        dNr2(2,3:3:r3,6) = -(1+xsi)/8;
        dNr2(2,3:3:r3,7) = (1+xsi)/8;
        dNr2(2,3:3:r3,8) = (1-xsi)/8;
        
        dNr2(3,1:3:r3,1) = (1-eta)/8;
        dNr2(3,1:3:r3,2) = -(1-eta)/8;
        dNr2(3,1:3:r3,3) = -(1+eta)/8;
        dNr2(3,1:3:r3,4) = (1+eta)/8;
        dNr2(3,1:3:r3,5) = -(1-eta)/8;
        dNr2(3,1:3:r3,6) = (1-eta)/8;
        dNr2(3,1:3:r3,7) = (1+eta)/8;
        dNr2(3,1:3:r3,8) = -(1+eta)/8;
        
        dNr2(3,2:3:r3,1) = (1-xsi)/8;
        dNr2(3,2:3:r3,2) = (1+xsi)/8;
        dNr2(3,2:3:r3,3) = -(1+xsi)/8;
        dNr2(3,2:3:r3,4) = -(1-xsi)/8;
        dNr2(3,2:3:r3,5) = -(1-xsi)/8;
        dNr2(3,2:3:r3,6) = -(1+xsi)/8;
        dNr2(3,2:3:r3,7) = (1+xsi)/8;
        dNr2(3,2:3:r3,8) = (1-xsi)/8;
        
        dNr2(3,3:3:r3,1) = 0;
        dNr2(3,3:3:r3,2) = 0;
        dNr2(3,3:3:r3,3) = 0;
        dNr2(3,3:3:r3,4) = 0;
        dNr2(3,3:3:r3,5) = 0;
        dNr2(3,3:3:r3,6) = 0;
        dNr2(3,3:3:r3,7) = 0;
        dNr2(3,3:3:r3,8) = 0;
    else
        error('not implemented');
    end
else
    error('dimension not implemented')    
end

%% output
out.gaussWeight = gaussWeight;
out.gaussPoint = gaussPoint;
out.N = N;
out.dNr = dNr;
out.dNr2 = dNr2;

end

