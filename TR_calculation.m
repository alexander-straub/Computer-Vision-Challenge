function [T, R, lambda, P1] = TR_reconstruction(T1, T2, R1, R2, correspondences, K)
    %% Preparation
    % calculate N
    dim = size(correspondences);
    N = dim(2)
    
    % create cells
    T_cell = {T1,T2,T1,T2};
    R_cell = {R1,R1,R2,R2};
    d_cell = {zeros(N,2),zeros(N,2),zeros(N,2),zeros(N,2)};
    
    % construct homogenous uncalibrated coordinates
    x1_uc = [correspondences(1:2,:);ones(1,N)];
    x2_uc = [correspondences(3:4,:);ones(1,N)];
    
    % calibrate homogenous coordinates
    x1 = inv(K)*x1_uc;
    x2 = inv(K)*x2_uc;
    
    %% Reconstruction

    dim = size(x1);
    for i = 1:4
        R = R_cell{i};
        T = T_cell{i};
        
        [M1,M2] = build_matrices(x1,x2,T,R,dim);
        
        [U1,S1,V1] = svd(M1);
        [U2,S2,V2] = svd(M2);
        
        d1 = V1(:,end);
        d2 = V2(:,end);
        
        d1 = d1/d1(end);
        d2 = d2/d2(end);

        d_cell{i} = [d1(1:end-1),d2(1:end-1)];
    end
    
    val1 = find(d_cell{1} > 0);
    val2 = find(d_cell{2} > 0);
    val3 = find(d_cell{3} > 0);
    val4 = find(d_cell{4} > 0);
    [M,I] = max([length(val1),length(val2),length(val3),length(val4)]);
    
    T = T_cell{I};
    R = R_cell{I};
    lambda = d_cell{I};
    
    [M1,M2] = build_matrices(x1,x2,T2,R2,dim);

    %% Calculation and visualization of the 3D points and the cameras
     dim = size(x1);
    %camC1 = [-0.2,0.2,0.2,-0.2;0.2,0.2,-0.2,-0.2;1,1,1,1];
    %camC2 = R'*camC1 - repmat(R'*T,1,4);
    idx = [1:dim(2)];
    P1 = x1.*lambda(:,1)';
    figure 
    hold on
    scatter3(P1(1,:),P1(2,:),P1(3,:),'o');
    text(P1(1,:),P1(2,:),P1(3,:),string(idx));
    %plot3([camC1(1,:),camC1(1,1)],[camC1(2,:),camC1(2,1)],[camC1(3,:),camC1(3,1)],'blue');
    %plot3([camC2(1,:),camC2(1,1)],[camC2(2,:),camC2(2,1)],[camC2(3,:),camC2(3,1)],'red');
    %text(camC1(1,4),camC1(2,4),camC1(3,4),'Cam1','Color','blue');
    %text(camC2(1,4),camC2(2,4),camC2(3,4),'Cam2','Color','red');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');  
    grid on;
    campos([43,-22,-87]);
    camup([0,-1,0]);
    
    
        
end

function [M1,M2] = build_matrices(x1, x2, T, R, dim)
        
        M1_diag = cross(x2,R*x1);
        M2_diag = cross(x1,R'*x2);
        
        %
        T_mat = repmat(T,1,dim(2));
        x2_T = cross(x2,T_mat);
        x2_T_s = reshape(x2_T,[],1);
        
        %
        R_T_mat = R'*T_mat;        
        x1_RT = -cross(x1, R_T_mat);
        x1_RT_s = reshape(x1_RT,[],1);
        
        % initialize M
        M1 = zeros(3*dim(2),dim(2)+1);
        M2 = zeros(3*dim(2),dim(2)+1);
        
        for k = 1:dim(2)
            M1((3*k-2):(3*k),k) = M1_diag(:,k);
            M2((3*k-2):(3*k),k) = M2_diag(:,k);
        end
         
        M1(:,end) = x2_T_s;
        M2(:,end) = x1_RT_s;
end