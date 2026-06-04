function [v] = three_point_system(E,T,x1,x2)
    
    x2_hat = hat(x2);
    T_hat = hat(T);
    b = x2_hat*T_hat'*E*x1;
    A = repmat(x2_hat*T,1,3).* repmat(x1',3,1);
    v = A\b;

end