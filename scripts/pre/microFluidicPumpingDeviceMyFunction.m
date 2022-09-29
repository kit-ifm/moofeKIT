function out = microFluidicPumpingDeviceMyFunction(t)
phi0 = 7.5e8;
out = 0;
% C = {@sin,@sqrt,@pow2};
% for k = 1:numel(C)
%     out = out + C{k}(x);
% end
timeFunction2 = {};
for k=1:4
    tk0 = (k-1)*32;
    tk1 = tk0 + 4;
    tk2 = tk1 + 4;
    tk3 = tk2 + 4;
    tk4 = tk3 + 4;
    timeFunction2{k} = @(t) (phi0*sin(pi/2*(t-tk0)/(tk1-tk0))).*(t >= tk0).*(t <= tk1) + ...
                            phi0*(t > tk1).*(t <= tk2) + ...
                            (phi0*cos(pi/2*(t-tk2)/(tk3-tk2))).*(t > tk2).*(t <= tk3);
     out = out + timeFunction2{k}(t);
end
end