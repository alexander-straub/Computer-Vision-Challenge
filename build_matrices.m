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
