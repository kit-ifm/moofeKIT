function B = BMatrix(dNX,varargin)
% nodal operator matrix for linear (empty varargin) and nonlinear (varargin
% = F) case
ii = 1;
while ii <= numel(varargin)
    if strcmpi(varargin{ii},'mapType')
        mapType = varargin{ii+1};
        ii = ii + 1;
    elseif strcmpi(varargin{ii},'mapVoigtObject')
        mapVoigtObject = varargin{ii+1};
        mapType = mapVoigtObject.mapType;
        ii = ii + 1;
    else
        F = varargin{1};
    end
    ii = ii + 1;
end
dimension = size(dNX,1);
numberOfNodes = size(dNX,2);
numberOfDisplacementDofs = numberOfNodes*dimension;
numberOfSymmetricVoigtDofs = dimension + (dimension^2-dimension)/2;
B = zeros(numberOfSymmetricVoigtDofs,numberOfDisplacementDofs);
if ~(exist('F', 'var') == 1)
    if ~exist('mapType', 'var')
        mapType = 'symmetric';
    end
    if strcmpi(mapType,'symmetric')
        if dimension == 1
            B = dNX;
        elseif dimension == 2
            B(1,1:dimension:end) = dNX(1,:);
            B(2,2:dimension:end) = dNX(2,:);
            B(3,1:dimension:end) = dNX(2,:);
            B(3,2:dimension:end) = dNX(1,:);
        elseif dimension == 3
            B(1,1:dimension:end) = dNX(1,:);
            B(2,2:dimension:end) = dNX(2,:);
            B(3,3:dimension:end) = dNX(3,:);
            B(4,1:dimension:end) = dNX(2,:);
            B(4,2:dimension:end) = dNX(1,:);
            B(5,2:dimension:end) = dNX(3,:);
            B(5,3:dimension:end) = dNX(2,:);
            B(6,1:dimension:end) = dNX(3,:);
            B(6,3:dimension:end) = dNX(1,:);
        end
    elseif strcmpi(mapType,'unsymmetric')
        B = zeros(dimension^2,size(dNX,2)*dimension);
        if dimension == 2
            % B(1,1:dimension:end) = dNX(1,:);
            % B(2,2:dimension:end) = dNX(2,:);
            % B(3,1:dimension:end) = dNX(2,:);
            % B(4,2:dimension:end) = dNX(1,:);
            % FIXME discuss!
            B(1,1:dimension:end) = dNX(1,:);
            B(2,2:dimension:end) = dNX(2,:);
            B(3,1:dimension:end) = 1/2*dNX(2,:);
            B(3,2:dimension:end) = 1/2*dNX(1,:);
            B(4,1:dimension:end) = 1/2*dNX(2,:);
            B(4,2:dimension:end) = 1/2*dNX(1,:);
        elseif dimension == 3
            % B(1,1:dimension:end) = dNX(1,:);
            % B(2,2:dimension:end) = dNX(2,:);
            % B(3,3:dimension:end) = dNX(3,:);
            % B(4,1:dimension:end) = dNX(2,:);
            % B(5,2:dimension:end) = dNX(3,:);
            % B(6,1:dimension:end) = dNX(3,:);
            % B(7,2:dimension:end) = dNX(1,:);
            % B(8,3:dimension:end) = dNX(2,:);
            % B(9,3:dimension:end) = dNX(1,:);
            % FIXME discuss!
            B(1,1:dimension:end) = dNX(1,:);
            B(2,2:dimension:end) = dNX(2,:);
            B(3,3:dimension:end) = dNX(3,:);
            B(4,1:dimension:end) = 1/2*dNX(2,:);
            B(4,2:dimension:end) = 1/2*dNX(1,:);
            B(5,2:dimension:end) = 1/2*dNX(3,:);
            B(5,3:dimension:end) = 1/2*dNX(2,:);
            B(6,1:dimension:end) = 1/2*dNX(3,:);
            B(6,3:dimension:end) = 1/2*dNX(1,:);
            B(7,1:dimension:end) = 1/2*dNX(2,:);
            B(7,2:dimension:end) = 1/2*dNX(1,:);
            B(8,2:dimension:end) = 1/2*dNX(3,:);
            B(8,3:dimension:end) = 1/2*dNX(2,:);
            B(9,1:dimension:end) = 1/2*dNX(3,:);
            B(9,3:dimension:end) = 1/2*dNX(1,:);
        end
    end
else
    if ~exist('mapType', 'var')
        mapType = 'symmetric';
    end
    if strcmpi(mapType,'symmetric')
        if dimension == 2
            B(1,1:dimension:end) = F(1,1)*dNX(1,:);
            B(1,2:dimension:end) = F(2,1)*dNX(1,:);
            B(2,1:dimension:end) = F(1,2)*dNX(2,:);
            B(2,2:dimension:end) = F(2,2)*dNX(2,:);
            B(3,1:dimension:end) = F(1,1)*dNX(2,:) + F(1,2)*dNX(1,:);
            B(3,2:dimension:end) = F(2,1)*dNX(2,:) + F(2,2)*dNX(1,:);
        elseif dimension == 3
            B(1,1:dimension:end) = F(1,1)*dNX(1,:);
            B(1,2:dimension:end) = F(2,1)*dNX(1,:);
            B(1,3:dimension:end) = F(3,1)*dNX(1,:);
            B(2,1:dimension:end) = F(1,2)*dNX(2,:);
            B(2,2:dimension:end) = F(2,2)*dNX(2,:);
            B(2,3:dimension:end) = F(3,2)*dNX(2,:);
            B(3,1:dimension:end) = F(1,3)*dNX(3,:);
            B(3,2:dimension:end) = F(2,3)*dNX(3,:);
            B(3,3:dimension:end) = F(3,3)*dNX(3,:);
            B(4,1:dimension:end) = F(1,1)*dNX(2,:) + F(1,2)*dNX(1,:);
            B(4,2:dimension:end) = F(2,1)*dNX(2,:) + F(2,2)*dNX(1,:);
            B(4,3:dimension:end) = F(3,1)*dNX(2,:) + F(3,2)*dNX(1,:);
            B(5,1:dimension:end) = F(1,2)*dNX(3,:) + F(1,3)*dNX(2,:);
            B(5,2:dimension:end) = F(2,2)*dNX(3,:) + F(2,3)*dNX(2,:);
            B(5,3:dimension:end) = F(3,2)*dNX(3,:) + F(3,3)*dNX(2,:);
            B(6,1:dimension:end) = F(1,1)*dNX(3,:) + F(1,3)*dNX(1,:);
            B(6,2:dimension:end) = F(2,1)*dNX(3,:) + F(2,3)*dNX(1,:);
            B(6,3:dimension:end) = F(3,1)*dNX(3,:) + F(3,3)*dNX(1,:);
        end
    elseif strcmpi(mapType,'unsymmetric')
        if dimension == 2
            B(1,1:dimension:end) = F(1,1)*dNX(1,:);
            B(1,2:dimension:end) = F(2,1)*dNX(1,:);
            B(2,1:dimension:end) = F(1,2)*dNX(2,:);
            B(2,2:dimension:end) = F(2,2)*dNX(2,:);
            B(3,1:dimension:end) = 1/2*(F(1,1)*dNX(2,:) + F(1,2)*dNX(1,:));
            B(3,2:dimension:end) = 1/2*(F(2,1)*dNX(2,:) + F(2,2)*dNX(1,:));
            B(4,1:dimension:end) = 1/2*(F(1,1)*dNX(2,:) + F(1,2)*dNX(1,:));
            B(4,2:dimension:end) = 1/2*(F(2,1)*dNX(2,:) + F(2,2)*dNX(1,:));
        elseif dimension == 3
            B(1,1:dimension:end) = F(1,1)*dNX(1,:);
            B(1,2:dimension:end) = F(2,1)*dNX(1,:);
            B(1,3:dimension:end) = F(3,1)*dNX(1,:);
            B(2,1:dimension:end) = F(1,2)*dNX(2,:);
            B(2,2:dimension:end) = F(2,2)*dNX(2,:);
            B(2,3:dimension:end) = F(3,2)*dNX(2,:);
            B(3,1:dimension:end) = F(1,3)*dNX(3,:);
            B(3,2:dimension:end) = F(2,3)*dNX(3,:);
            B(3,3:dimension:end) = F(3,3)*dNX(3,:);
            B(4,1:dimension:end) = 1/2*(F(1,1)*dNX(2,:) + F(1,2)*dNX(1,:));
            B(4,2:dimension:end) = 1/2*(F(2,1)*dNX(2,:) + F(2,2)*dNX(1,:));
            B(4,3:dimension:end) = 1/2*(F(3,1)*dNX(2,:) + F(3,2)*dNX(1,:));
            B(5,1:dimension:end) = 1/2*(F(1,2)*dNX(3,:) + F(1,3)*dNX(2,:));
            B(5,2:dimension:end) = 1/2*(F(2,2)*dNX(3,:) + F(2,3)*dNX(2,:));
            B(5,3:dimension:end) = 1/2*(F(3,2)*dNX(3,:) + F(3,3)*dNX(2,:));
            B(6,1:dimension:end) = 1/2*(F(1,1)*dNX(3,:) + F(1,3)*dNX(1,:));
            B(6,2:dimension:end) = 1/2*(F(2,1)*dNX(3,:) + F(2,3)*dNX(1,:));
            B(6,3:dimension:end) = 1/2*(F(3,1)*dNX(3,:) + F(3,3)*dNX(1,:));
            B(7,1:dimension:end) = 1/2*(F(1,1)*dNX(2,:) + F(1,2)*dNX(1,:));
            B(7,2:dimension:end) = 1/2*(F(2,1)*dNX(2,:) + F(2,2)*dNX(1,:));
            B(7,3:dimension:end) = 1/2*(F(3,1)*dNX(2,:) + F(3,2)*dNX(1,:));
            B(8,1:dimension:end) = 1/2*(F(1,2)*dNX(3,:) + F(1,3)*dNX(2,:));
            B(8,2:dimension:end) = 1/2*(F(2,2)*dNX(3,:) + F(2,3)*dNX(2,:));
            B(8,3:dimension:end) = 1/2*(F(3,2)*dNX(3,:) + F(3,3)*dNX(2,:));
            B(9,1:dimension:end) = 1/2*(F(1,1)*dNX(3,:) + F(1,3)*dNX(1,:));
            B(9,2:dimension:end) = 1/2*(F(2,1)*dNX(3,:) + F(2,3)*dNX(1,:));
            B(9,3:dimension:end) = 1/2*(F(3,1)*dNX(3,:) + F(3,3)*dNX(1,:));
        end
    end
end