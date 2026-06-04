function Cake = cake(min_dist)
    % The cake function creates a "cake matrix" that contains a circular set-up of zeros
    % and fills the rest of the matrix with ones. 
    % This function can be used to eliminate all potential features around a stronger feature
    % that don't meet the minimal distance to this respective feature.
    x = 1:min_dist*2+1;
    y = 1:min_dist*2+1;
    [xx yy] = meshgrid(x,y);
    Cake = ones(size(xx));
    Cake(((xx-(min_dist+1)).^2+(yy-(min_dist+1)).^2)<=(min_dist)^2)=0;
    Cake = logical(Cake);
end