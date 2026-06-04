function [feat_corn,feat_edge,feat_homo] = harris_detector(input_image, varargin)
    % This function implements a Harris detector that extracts features
    % from the input_image.

    %% Input parser
    % create default values
    default_segment_length = 15;
    default_k = 0.05;
    default_tau = 10^6;
    default_do_plot = false(1);
    default_min_dist = 20;
    default_tile_size = 200;
    default_N = 5;

    % create parser
    p = inputParser;

    % create validation condition
    valid_segment_length = @(x) isnumeric(x) && (mod(x,2)~=0) && (x > 1);
    valid_k = @(x) isnumeric(x) && (x >= 0) && (x <= 1);
    valid_tau = @(x) isnumeric(x) && (x > 0);
    valid_do_plot = @(x) islogical(x);
    valid_input_image = @(x) ismatrix(x);
    valid_min_dist = @(x) isnumeric(x) && (x >= 1);
    valid_tile_size = @(x) isnumeric(x) || ALL(isvector(x), size(x) == [1 2]);
    valid_N = @(x) isnumeric(x) && (x >= 1);
   
    % check parameter
    addRequired(p,'input_image',valid_input_image);
    addParameter(p,'segment_length',default_segment_length,valid_segment_length);
    addParameter(p,'k',default_k,valid_k);
    addParameter(p,'tau',default_tau,valid_tau);
    addParameter(p,'do_plot',default_do_plot,valid_do_plot);
    addParameter(p,'min_dist',default_min_dist,valid_min_dist);
    addParameter(p,'tile_size',default_tile_size,valid_tile_size);
    addParameter(p,'N',default_N,valid_N);

    % parse arguments
    parse(p,input_image,varargin{:});

    % assign variables
    segment_length = p.Results.segment_length;
    k = p.Results.k;
    tau = p.Results.tau;
    do_plot = p.Results.do_plot;
    min_dist = p.Results.min_dist;
    N = p.Results.N;
    
    % check for input type of tile_size and assign accordingly
    if isscalar(p.Results.tile_size) 
        tile_size = [p.Results.tile_size, p.Results.tile_size];
    else 
        tile_size = p.Results.tile_size;
    end

     % Check if input_image is grayscale image
    dim = size(input_image);
    if length(dim) ~= 2
        error("Image format has to be NxMx1");
    end
    
    % Approximation of the image gradient
    gray = double(input_image);
    [Ix, Iy] = sobel_xy(gray);
    
    %% Weighting
    
    % parameter for gauss-function
    std = segment_length;
    mean = (segment_length+1)/2;
    steps = [1:segment_length];
    
    % gauss-vector of size segment_length
    gauss = (1/(sqrt(2*pi)*std))*exp((-1/2)*((steps-mean)/std).^2);

    % normalize gauss vector to sum up to 1
    w = gauss./sum(gauss);
    
    % calculate weight-matrix W
    W = w' * w;
    
    % Harris Matrix G
    % calculate gradient matrices
    grad_xx = Ix .* Ix;
    grad_xy = Ix .* Iy;
    grad_yy = Iy .* Iy;
    
    % calculate Harris Matrix G
    G11_uncut = conv2(grad_xx, W);
    G12_uncut = conv2(grad_xy, W);
    G22_uncut = conv2(grad_yy, W);

    
    % cutoff unwanted border from convolution
    cutoff = (segment_length+1)/2;
    G11 = G11_uncut(cutoff:dim(1)+cutoff-1, cutoff:dim(2)+cutoff-1);
    G12 = G12_uncut(cutoff:dim(1)+cutoff-1, cutoff:dim(2)+cutoff-1);
    G22 = G22_uncut(cutoff:dim(1)+cutoff-1, cutoff:dim(2)+cutoff-1);

    %% Feature extraction with the Harris measurement
    % initialize variables
    dim = size(G11);
    H = zeros(dim(1),dim(2));
    corners = zeros(dim(1),dim(2));
    %edges = zeros(dim(1),dim(2));
    %homogene = zeros(dim(1),dim(2));
    
    % calculate Harris-matrice
    for k = 1:dim(1)
        for l = 1:dim(2)
            G_x = [G11(k,l) G12(k,l); G12(k,l) G22(k,l)];
            H(k,l) = det(G_x) - 0.05 *  (trace(G_x))^2;
        end
    end
    
    %H_r = round(H,2,'significant');
    %histogram(H_r)

    %tau1 = 0.75*10^8;
    %tau2 = -0.075*10^8;
    %tau3 = 10^4;

    % create corner-matrice with elements of H > threshold
    corners(H > tau) = H (H > tau);
    % create corner-matrice with elements of H > threshold
    %edges(H < tau2) = H (H < tau2);
    % create corner-matrice with elements of H > threshold
    %homogene(abs(H) < tau3) = H (abs(H) < tau3);
    
    % cutoff pixels at the border (high values at the border due to calculation method of gradient)
    cutoff = ceil(segment_length/2); 
    corners(1:cutoff,:) = 0;
    corners(:,1:cutoff) = 0;
    corners(dim(1)-cutoff:dim(1),:) = 0;
    corners(:,dim(2)-cutoff:dim(2)) = 0;
    % edges(1:cutoff,:) = 0;
    % edges(:,1:cutoff) = 0;
    % edges(dim(1)-cutoff:dim(1),:) = 0;
    % edges(:,dim(2)-cutoff:dim(2)) = 0;
    % homogene(1:cutoff,:) = 0;
    % homogene(:,1:cutoff) = 0;
    % homogene(dim(1)-cutoff:dim(1),:) = 0;
    % homogene(:,dim(2)-cutoff:dim(2)) = 0;
   
    % create feature matrice with coordinates of pixels with high feature values (corners)
    [row, col] = find(corners);
    % [row1, col1] = find(edges);
    % [row2, col2] = find(homogene);
    feat_corner = [col'; row'];
    % feat_edges = [col1'; row1'];
    % feat_homog = [col2'; row2'];

    %% Plot
    if do_plot 
        figure;
        % subplot(2,2,1);
        imshow(uint8(input_image));
        hold on
        plot(feat_corner(1,:), feat_corner(2,:),"o");
        hold off
        % subplot(2,2,2);
        % imshow(uint8(input_image));
        % hold on
        % plot(feat_edges(1,:), feat_edges(2,:),"o");
        % hold off
        % subplot(2,2,3);
        % imshow(uint8(input_image));
        % hold on
        % plot(feat_homog(1,:), feat_homog(2,:),"o");
        % hold off
    end

    %% Feature preparation
    dim_corners = size(corners);
    % initialize zero-matrice for border
    border_corners = zeros(dim_corners(1)+2*min_dist, dim_corners(2)+2*min_dist);
    % border_edges = border_corners;
    % border_homog = border_corners;
    
    % put values of corners into middle of border_corners;
    border_corners(min_dist+1:end-min_dist, min_dist+1:end-min_dist) = corners;
    % border_edges(min_dist+1:end-min_dist, min_dist+1:end-min_dist) = edges;
    % border_homog(min_dist+1:end-min_dist, min_dist+1:end-min_dist) = homogene;
    corners = border_corners;
    % edges = border_edges;
    % homogene = border_homog;
    
    % find value and coordinates of non-zero elements in corners
    [y, x, val] = find(corners);
    % [y1, x1, val1] = find(edges);
    % [y2, x2, val2] = find(homogene);
    
    % put values and coordinates in one matrice and sort it by decreasing values
    feat_corner_unsorted = [y, x, val];
    % feat_edges_unsorted = [y1, x1, val1];
    % feat_homog_unsorted = [y2, x2, val2];
    feat_corner_sorted = sortrows(feat_corner_unsorted,3,'descend');
    % feat_edges_sorted = sortrows(feat_edges_unsorted,3,'ascend');
    % feat_homog_sorted = sortrows(abs(feat_homog_unsorted),3,'ascend');
    
    % put index of sorted values in a column vector by calculating: x + x_dim * y
    %dim_border_corners = size(border_corners);
    %sorted_index = [features_sorted(:,2) + dim_border_corners(1) .* (features_sorted(:,1)-1) ];

    %% Accumulator array
    % calculate the number of tiles fitting into input_image
    dim_ii = size(input_image);
    N_x_tiles = ceil(dim_ii(2)/tile_size(2));
    N_y_tiles = ceil(dim_ii(1)/tile_size(1));
    acc_array = zeros(N_y_tiles, N_x_tiles);
    
    % calculate the maximum number of features possible and initialize features
    total_N_feat_corn = min([length(val), N_y_tiles*N_x_tiles*N]);
    % total_N_feat_edge = min([length(val1), N_y_tiles*N_x_tiles*N]);
    % total_N_feat_homo = min([length(val2), N_y_tiles*N_x_tiles*N]);
    feat_corn = zeros(2, total_N_feat_corn);
    % feat_edge = zeros(2, total_N_feat_edge);
    % feat_homo = zeros(2, total_N_feat_homo);

    %% Feature detection with minimal distance and maximal number of features per tile
    % calculate coordinates of all features in sorted_index
    %dim_c = size(corners);
    %x_vec = ceil(sorted_index./dim_c(1));
    %coordinates = [x_vec, sorted_index-(x_vec-1).*dim_c(1)];
    
    [feat_corn, N_curr_f_c] = feat_per_tile(feat_corner_sorted,min_dist,tile_size,N,acc_array,feat_corn,corners);
    % [feat_edge, N_curr_f_e] = feat_per_tile(feat_edges_sorted,min_dist,tile_size,N,acc_array,feat_edge,edges);
    % [feat_homo, N_curr_f_h] = feat_per_tile(feat_homog_sorted,min_dist,tile_size,N,acc_array,feat_homo,homogene);
   
    % cutoff part of feature-vector containing only zeros
    feat_corn = feat_corn(1:2, 1:N_curr_f_c);
    % feat_edge = feat_edge(1:2, 1:N_curr_f_e);
    % feat_homo = feat_homo(1:2, 1:N_curr_f_h);
    
    % Plot Routine
    figure;
    imshow(input_image);
    hold on
    plot(feat_corn(1,:),feat_corn(2,:),'*blue');
    % plot(feat_edge(1,:),feat_edge(2,:),'*red');
    % plot(feat_homo(1,:),feat_homo(2,:),'*yellow');
    hold off
    feat_edge = 1;
    feat_homo = 1;

end
