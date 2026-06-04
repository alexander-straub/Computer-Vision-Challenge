function W = hat(w)
% This function implements the ^-operator.
% It converts a 3x1-vector into a skew symmetric matrix.
W = [0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0];
end
