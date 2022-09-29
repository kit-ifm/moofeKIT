%% solver
% parpool(24)
% 
warning off;
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:singularMatrix')
parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix')
% try
    dofObject = runNewton(setupObject,dofObject);
% catch
%     iiVector = [iiVector;ii];
% end
% iiVector
% end