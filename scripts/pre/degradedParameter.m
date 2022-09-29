function [pN1] = degradedParameter(p0,pInf,eRef,eps0,r,dissWorkN1)
 % Damage model: pN1 is chosen material constant from Poynting-Thomson-model
 %               p0, pInf are initial, final values of pN1
 %               eps0 = amplitude of strain
 %               eref is chosen reference e-modulus

 
 
 
 pN1=p0*exp(-r.*dissWorkN1/(eRef*eps0*eps0))+pInf;
end
