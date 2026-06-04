function [H] = four_point_alg_ess(correspondences, K)
    
    dim = size(correspondences);
    x1 = correspondences(4:6,:);
    x2 = correspondences(1:3,:);

    A = zeros(dim(2)*3,9);
    for i = 1:dim(2);
        x_2_hat = hat(x2(:,i));
        B = kron(x1(:,i),x_2_hat);
        A(3*i-2:3*i,:) = B';
    end

    [U_A,sigma_A,V_A] = svd(A);
    H_L = reshape(V_A(:,end),[3,3]);
    [U_H_L,sigma_H_L,V_H_L] = svd(H_L);
    sigma_H_L(2,2)
    H = (1/sigma_H_L(2,2))*H_L;

end