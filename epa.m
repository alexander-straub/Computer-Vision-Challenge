function EF = epa(correspondences, K)
    % Depending on whether a calibrating matrix 'K' is given,
    % this function calculates either the essential or the fundamental matrix
    % with the eight-point algorithm.

    N_corr = size(correspondences); 
    switch nargin
        case 2
            x1_unc = [correspondences(1:2,1:end); ones(1,N_corr(2))];
            x2_unc = [correspondences(3:4,1:end); ones(1,N_corr(2))];
            x1 = inv(K) * x1_unc;
            x2 = inv(K) * x2_unc;
        case 1
            x1 = [correspondences(1:2,1:end); ones(1,N_corr(2))];
            x2 = [correspondences(3:4,1:end); ones(1,N_corr(2))];
    end
    
    A = zeros(N_corr(2),9);
    for i = 1:N_corr(2)
        A(i,:) = kron(x1(:,i),x2(:,i))';
    end
    
    [U,S,V] = svd(A);

    %% Estimation of the matrices
    G_s = reshape(V(:,end),[3,3]);
    [Ug,Sg,Vg] = svd(G_s);
    switch nargin
        case 2
            Sg(1,1) = 1;
            Sg(2,2) = 1;
            Sg(3,3) = 0;
        case 1
            Sg(3,3) = 0;
    end
    EF = Ug * Sg * Vg';
    
end