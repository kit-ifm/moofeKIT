function out = selectStress(stressTensor,setupObject,dimension)
% This function selects the components out of the stress tensor to be
% plotted. It works also for nonsymetric stress tensors.
%
% 11.11.2016 F.Streich
component = setupObject.plotObject.stress.component;
type = setupObject.plotObject.stress.type;
stress = stressTensor.(type);
if dimension == 1
    if isempty(component)
        out = stress;
    elseif component == 1
        out = stress(1);
    elseif component == 2
        out = stress(2);
    elseif component == 3
        out = stress(3);
    else
        error("please select valid stress component for postprocessing!")
    end
elseif dimension == 2
    assert(ismember(component, [-1,11,12,21,22]), 'Third components are only in 3D availible')
    if component == -1 % von Mises stress
        out = real(sqrt(stress(1,1)^2 + stress(2,2)^2 - stress(1,1)*stress(2,2) + 3*stress(1,2)*stress(2,1)));
    elseif component == 11
        out = stress(1,1);
    elseif component == 22
        out = stress(2,2);
    elseif component == 12
        out = stress(1,2);
    elseif component == 21
        out = stress(2,1);
    end
elseif dimension == 3
    if component == -1 % von Mises stress
        out = real(sqrt(stress(1,1)^2 + stress(2,2)^2 + stress(3,3)^2 - ...
            stress(1,1)*stress(2,2) - stress(2,2)*stress(3,3) - stress(3,3)*stress(1,1) +...
            3*stress(1,2)*stress(2,1) + 3*stress(2,3)*stress(3,2) + 3*stress(1,3)*stress(3,1)));
    elseif component == 11
        out = stress(1,1);
    elseif component == 22
        out = stress(2,2);
    elseif component == 33
        out = stress(3,3);
    elseif component == 12
        out = stress(1,2);
    elseif component == 21
        out = stress(2,1);
    elseif component == 13
        out = stress(1,3);
    elseif component == 31
        out = stress(3,1);
    elseif component == 23
        out = stress(2,3);
    elseif component == 32
        out = stress(3,2);
    elseif strcmpi(component, 'maxPrincipalStress')
        out = max(eig(stress));
    elseif strcmpi(component, 'minPrincipalStress')
        out = min(eig(stress));
    end
end
end