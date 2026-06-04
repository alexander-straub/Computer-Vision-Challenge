function [correspondences_robust, largest_set_F] = F_ransac(correspondences, varargin)
    % This function implements the RANSAC algorithm to determine 
    % robust corresponding image points

    default_epsilon = 0.5;
    default_p = 0.5;
    default_tolerance = 0.01;
    % create parser
    parser = inputParser;

    % create validation condition
    valid_epsilon = @(x) isnumeric(x) && (x > 0) && (x < 1);
    valid_p = @(x) isnumeric(x) &&  (x > 0) && (x < 1);
    valid_tolerance = @(x) isnumeric(x);
   
    % check parameter
    addRequired(parser,'correspondences');
    addParameter(parser,'epsilon',default_epsilon,valid_epsilon);
    addParameter(parser,'p',default_p,valid_p);
    addParameter(parser,'tolerance',default_tolerance,valid_tolerance);

    % parse arguments
    parse(parser,correspondences,varargin{:});

    % assign variables
    correspondences = parser.Results.correspondences;
    epsilon = parser.Results.epsilon;
    p = parser.Results.p;
    tolerance = parser.Results.tolerance;
    
    dim = size(correspondences);
    x1_pixel = [correspondences(1:2,1:end); ones(1,dim(2))];
    x2_pixel = [correspondences(3:4,1:end); ones(1,dim(2))];

    %% RANSAC algorithm preparation
    largest_set_size = 0;
    largest_set_dist = Inf;
    largest_set_F = zeros(3);
    k = 8;
    s = log(1-p)/log(1-(1-epsilon)^k);

     %% RANSAC algorithm
    n_corr = size(x1_pixel);
    for i = 1:s
        idx = ceil(rand(1,k)*n_corr(2));
        k_rand_corr = correspondences(:,idx);
        F = epa(k_rand_corr);
        sd = sampson_dist(F,x1_pixel,x2_pixel);
        [row,col,curr_set] = find(sd < tolerance);
        curr_set_size = length(curr_set);
        curr_set_dist = sum(curr_set,'all');
        if curr_set_size > largest_set_size
            largest_set_size = curr_set_size;
            largest_set_dist = curr_set_dist;
            largest_set_F = F;
            idx_largest_set = col;
        elseif curr_set_size == largest_set_size
            if curr_set_dist < largest_set_dist
                largest_set_size = curr_set_size;
                largest_set_dist = curr_set_dist;
                largest_set_F = F;
                idx_largest_set = col;
            end
        end
    end
    
    %correspondences_robust = [x1_pixel(1:2,idx_largest_set); x2_pixel(1:2,idx_largest_set)];
    correspondences_robust = correspondences(:,idx_largest_set);
end