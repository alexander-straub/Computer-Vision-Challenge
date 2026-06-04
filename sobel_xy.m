function [Fx, Fy] = sobel_xy(input_image)
% This function uses the Sobel filter to calculate
% the image gradient of a grayscale image in x- and y- direction.
    S_x = [1 0 -1;2 0 -2;1 0 -1];
    S_y = [1 2 1;0 0 0;-1 -2 -1];
    Fx = conv2(input_image,S_x);
    Fy = conv2(input_image,S_y);
    Fx = Fx(2:end-1,2:end-1);
    Fy = Fy(2:end-1,2:end-1);      
end