function dN_X_I = computedN_X_I(dN_xi_k_I, J, k)
%COMPUTEDN_X_I Summary of this function goes here
%   Detailed explanation goes here

dNr = reshape(dN_xi_k_I(:, k, :), size(dN_xi_k_I, 1), size(dN_xi_k_I, 3));

    if size(J, 1) == size(J, 2)
        dN_X_I = J' \ dNr;
    elseif size(J, 2) == 1 && size(J, 1) ~= 1
        detJ = norm(J);
        dN_X_I = detJ^(-1)*dNr;
    else
        error('not implemented yet.')
    end

end
