function [features, N_curr_f] = feat_per_tile(feat_sorted,min_dist,tile_size,N,acc_array,features,featmatrix)
     
    % initialize variables
    N_curr_f = 0;                   % number of current chosen features
    dim_f = size(features);         % dimensions of feature vector
    max_f = dim_f(2);               % maximum possible amount of features 
    Cake = cake (min_dist);         % quadratic matrice of size min_dist*2+1 
    
    % iterate through all features in val
    for idx = 1:max_f
        
        % coordinates of current feature
        x = feat_sorted(idx,2);
        y = feat_sorted(idx,1);

        % check if feature got eliminated by Cake
        if featmatrix(y,x) ~= 0
            
            % determine, which tile the current feature is part of
            x_n = ceil((x-min_dist)/tile_size(2));
            y_n = ceil((y-min_dist)/tile_size(1));
            
            % check if more features are allowed in this tile
            if acc_array(y_n, x_n) < N
                
                % increase number of features in current tile by one 
                acc_array(y_n, x_n) = acc_array(y_n, x_n) + 1;
                
                % calculate the number of current features over all tiles combined
                N_curr_f = sum(acc_array,"all");
                
                % write coordinates into feature array 
                features(1, N_curr_f) = feat_sorted(idx,2);
                features(2, N_curr_f) = feat_sorted(idx,1);
                                                
                % eliminate surrounding features
                featmatrix(y-min_dist : y+min_dist, x-min_dist : x+min_dist) = ...
                featmatrix(y-min_dist : y+min_dist, x-min_dist : x+min_dist) .* Cake;
            end       
        end
        
        % stop if maximum possible amount of features were found 
        if N_curr_f == max_f
            disp('finished');
            break;
        end
    end
end