function [H] = four_point_alg(matchedpoints, matchedpoints2, K)
    dim_mp = size(matchedpoints.Location);
    matched_points1_z = [matchedpoints.Location'; ones(1,dim_mp(1))];
    matched_points2_z = [matchedpoints2.Location'; ones(1,dim_mp(1))];

    matched_points_cal = inv(K)*matched_points1_z;
    matched_points2_cal = inv(K)*matched_points2_z;
    A = zeros(dim_mp(1)*3,9);
    for i = 1:dim_mp(1)
        x_2_hat = hat(matched_points2_cal(:,i));
        B = [matched_points_cal(1,i)*x_2_hat; matched_points_cal(2,i)*x_2_hat; matched_points_cal(3,i)*x_2_hat];
        A(3*i-2:3*i,:) = B';
    end

    [U_A,sigma_A,V_A] = svd(A);
    H_L = reshape(V_A(:,end),[3,3]);
    [U_H_L,sigma_H_L,V_H_L] = svd(H_L);
    H = (1/sigma_H_L(2,2))*H_L;

end