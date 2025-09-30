function [C,CUnsymmetricVoigt,CSymmetricVoigt] = voigtTensorComputation(A,B)
    C = zeros(3,3,3,3);
    for I = 1:3
        for J = 1:3
            for K = 1:3
                for L = 1:3
                    % C(I,J,K,L) = 1/2*(A(I,K)*B(L,J)+A(I,L)*B(K,J));
                    % C(I,J,K,L) = A(I,K)*B(L,J)+A(I,L)*B(K,J);
                    % C(I,J,K,L) = kroneckerDeltaFunction(I,K)*kroneckerDeltaFunction(J,L) + kroneckerDeltaFunction(I,L)*kroneckerDeltaFunction(J,K);
                    % for P = 1:3
                    %     C(I,J,K,L) = C(I,J,K,L) + leviCevitaFunction(I,P,K)*leviCevitaFunction(J,P,L);
                    % end
                    for P = 1:3
                        for M = 1:3
                            C(I,J,K,L) = C(I,J,K,L) + leviCevitaFunction(I,M,K)*leviCevitaFunction(J,P,L)*kroneckerDeltaFunction(M,P);
                        end
                    end
                end
            end
        end
    end
    % TODO: build an operator for different Voigt matrix index notation instead of just one for the following fourth order Voigt tensor
    CUnsymmetricVoigt = [  C(1,1,1,1) C(1,1,2,2) C(1,1,3,3) C(1,1,1,2) C(1,1,2,3) C(1,1,1,3) C(1,1,2,1) C(1,1,3,2) C(1,1,3,1);
                C(2,2,1,1) C(2,2,2,2) C(2,2,3,3) C(2,2,1,2) C(2,2,2,3) C(2,2,1,3) C(2,2,2,1) C(2,2,3,2) C(2,2,3,1);
                C(3,3,1,1) C(3,3,2,2) C(3,3,3,3) C(3,3,1,2) C(3,3,2,3) C(3,3,1,3) C(3,3,2,1) C(3,3,3,2) C(3,3,3,1);
                C(1,2,1,1) C(1,2,2,2) C(1,2,3,3) C(1,2,1,2) C(1,2,2,3) C(1,2,1,3) C(1,2,2,1) C(1,2,3,2) C(1,2,3,1);
                C(2,3,1,1) C(2,3,2,2) C(2,3,3,3) C(2,3,1,2) C(2,3,2,3) C(2,3,1,3) C(2,3,2,1) C(2,3,3,2) C(2,3,3,1);
                C(1,3,1,1) C(1,3,2,2) C(1,3,3,3) C(1,3,1,2) C(1,3,2,3) C(1,3,1,3) C(1,3,2,1) C(1,3,3,2) C(1,3,3,1);
                C(2,1,1,1) C(2,1,2,2) C(2,1,3,3) C(2,1,1,2) C(2,1,2,3) C(2,1,1,3) C(2,1,2,1) C(2,1,3,2) C(2,1,3,1);
                C(3,2,1,1) C(3,2,2,2) C(3,2,3,3) C(3,2,1,2) C(3,2,2,3) C(3,2,1,3) C(3,2,2,1) C(3,2,3,2) C(3,2,3,1);
                C(3,1,1,1) C(3,1,2,2) C(3,1,3,3) C(3,1,1,2) C(3,1,2,3) C(3,1,1,3) C(3,1,2,1) C(3,1,3,2) C(3,1,3,1)];
    CSymmetricVoigt = [  C(1,1,1,1) C(1,1,2,2) C(1,1,3,3) C(1,1,1,2) C(1,1,2,3) C(1,1,1,3);
                C(2,2,1,1) C(2,2,2,2) C(2,2,3,3) C(2,2,1,2) C(2,2,2,3) C(2,2,1,3);
                C(3,3,1,1) C(3,3,2,2) C(3,3,3,3) C(3,3,1,2) C(3,3,2,3) C(3,3,1,3);
                C(1,2,1,1) C(1,2,2,2) C(1,2,3,3) C(1,2,1,2) C(1,2,2,3) C(1,2,1,3);
                C(2,3,1,1) C(2,3,2,2) C(2,3,3,3) C(2,3,1,2) C(2,3,2,3) C(2,3,1,3);
                C(1,3,1,1) C(1,3,2,2) C(1,3,3,3) C(1,3,1,2) C(1,3,2,3) C(1,3,1,3)];
end      

function delta = kroneckerDeltaFunction(I,J)
if I == J
    delta = 1;
else
    delta = 0;
end
end

function espislon = leviCevitaFunction(I,J,K)
if I == J || J == K || I == K
    espislon = 0;
elseif (I==1 && J==2 && K==3) || (I==2 && J==3 && K==1) || (I==3 && J==1 && K==2)
    espislon = 1;
elseif (I==1 && J==3 && K==2) || (I==3 && J==2 && K==1) || (I==2 && J==1 && K==3)
    espislon = -1;
end
end
