function [T1, R1, T2, R2, U, V] = TR_from_E(E)
    % This function calculates the possible values for T and R 
    % from the essential matrix
    Rz_1 = [0,1,0;-1,0,0;0,0,1];
    Rz_2 = [0,-1,0;1,0,0;0,0,1];
    [U,S,V] = svd(E);
    if det(U) < 0
        U = U*[1,0,0;0,1,0;0,0,-1];
    end
    if det(V) < 0
        V = V*[1,0,0;0,1,0;0,0,-1];
    end    
    R1 = U*Rz_1'*V';
    R2 = U*Rz_2'*V';
    det(R1)
    det(R2)
    T1_hat = U*Rz_1*S*U';
    T2_hat = U*Rz_2*S*U';
    T1 = [-T1_hat(2,3);T1_hat(1,3);-T1_hat(1,2);];
    T2 = [-T2_hat(2,3);T2_hat(1,3);-T2_hat(1,2);];
end