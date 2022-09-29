function [H] = wedge(A,B)  
%% Creates tensor cross product.
    %
    % Syntax
    % A,B vector or matrix
    %
    % [H] = wedge(FAkt,FAkt)      
    %
    % Description
    %
    % Tensor cross product.
    
    % 16.01.2015 C. Bilgen
   
    if size(A) == size(B) 
        if  A == B
            H = zeros(size(A));
            H(1,:) = 2*cross(A(2,:),A(3,:));
            H(2,:) = 2*cross(A(3,:),A(1,:));
            H(3,:) = 2*cross(A(1,:),A(2,:));
        else
            H = zeros(size(A));
            H(1,:) = cross(A(2,:),B(3,:)) + cross(B(2,:),A(3,:));
            H(2,:) = cross(B(3,:),A(1,:)) + cross(A(3,:),B(1,:));
            H(3,:) = cross(A(1,:),B(2,:)) + cross(B(1,:),A(2,:));
        end
    elseif isvector(A) || isvector(B)
        if isvector(A) == 1 && isvector(B) == 0
            H = zeros(size(B));
            H(:,1) = cross(A,B(:,1));
            H(:,2) = cross(A,B(:,2));
            H(:,3) = cross(A,B(:,3));
        elseif isvector(A) == 0 && isvector(B) == 1
            H = zeros(size(A));
            H(1,:) = cross(A(1,:),B);       
            H(2,:) = cross(A(2,:),B);
            H(3,:) = cross(A(3,:),B);
        else 
            error('matrix dimensions do not agree!')
        end
       
    else
        error('Matrix dimensions do not agree!')
    end
    
    
end