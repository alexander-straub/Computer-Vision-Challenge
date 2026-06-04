function cor = point_correspondence(I1, I2, Ftp1, Ftp2, varargin)
    % In this function you are going to compare the extracted features of a stereo recording
    % with NCC to determine corresponding image points.
    
    %% Input parser
    % create default values
    default_window_length = 25;
    default_min_corr = 0.95;
    default_do_plot = false(1);
    % create parser
    p = inputParser;
    
    % create validation condition
    valid_I1 = @(x) ismatrix(x);
    valid_I2 = @(x) ismatrix(x);
    valid_Ftp2 = @(x) ismatrix(x);
    valid_Ftp1 = @(x) ismatrix(x);
    valid_window_length = @(x) isnumeric(x) && (mod(x,2)==1) && (x > 1);
    valid_min_corr = @(x) isnumeric(x) && (x > 0) && (x < 1);
    valid_do_plot = @(x) islogical(x);
    
    % check parameter
    addRequired(p,'I1',valid_I1);
    addRequired(p,'I2',valid_I2);
    addRequired(p,'Ftp1',valid_Ftp1);
    addRequired(p,'Ftp2',valid_Ftp2);
    addParameter(p,'window_length',default_window_length,valid_window_length);
    addParameter(p,'min_corr',default_min_corr,valid_min_corr);
    addParameter(p,'do_plot',default_do_plot,valid_do_plot);
    
    % parse arguments
    parse(p,I1,I2,Ftp1,Ftp2,varargin{:});
    
    % assign variables
    window_length = p.Results.window_length;
    min_corr = p.Results.min_corr;
    do_plot = p.Results.do_plot;
    Im1 = double(p.Results.I1);
    Im2 = double(p.Results.I2);

    %% Feature preparation
    

    border = floor(window_length/2);
    dim_image1 = size(I1);
    dim1 = size(Ftp1);
    dim_image2 = size(I2);
    dim2 = size(Ftp2);

    for i = 1:dim1(2)
        coord = Ftp1(:,i);
        if coord(1) <= border || coord(1) > dim_image1(2)-border
            Ftp1(1:2,i) = zeros(2,1);
        elseif coord(2) <= border || coord(2) > dim_image1(1)-border
            Ftp1(1:2,i) = zeros(2,1);
        else
            Ftp1(1:2,i) = coord;
        end
    end

    for j = 1:dim2(2)
        coord2 = Ftp2(:,j);
        if coord2(1) <= border || coord2(1) > dim_image2(2)-border
            Ftp2(1:2,j) = zeros(2,1);
        elseif coord2(2) <= border || coord2(2) > dim_image2(1)-border
            Ftp2(1:2,j) = zeros(2,1);
        else
            Ftp2(1:2,j) = coord2;
        end
    end

       
    [~,col1] = find(Ftp1(1,:));
    [~,col2] = find(Ftp2(1,:));
    
    Ftp1 = Ftp1(:,col1);
    Ftp2 = Ftp2(:,col2);
    

    %% Normalization
    % convert to double for calculation
    I1 = double(I1);
    I2 = double(I2);
    
    % get parameter and initialize variables
    N = window_length^2;
    dist = floor(window_length/2);
    Mat_feat_1 = zeros(N,length(Ftp1(1,:)));
    Mat_feat_2 = zeros(N,length(Ftp2(1,:)));
    
    % feature window operation for features of I1
    for i= 1:length(Ftp1(1,:))
        
        % coordinates of current feature
        x = Ftp1(1,i);
        y = Ftp1(2,i);
        
        % window around current feature
        W = I1(y-dist:y+dist,x-dist:x+dist);
        
        % mean, std of window and the normalized Window W_n
        W_mean = (sum(W,"all")/N)*ones(window_length);
        W_std = sqrt((1/(N-1))*trace((W-W_mean)'*(W-W_mean)));
        W_n = (1/W_std)*(W-W_mean);
        
        % convert normalized Window into col-vector and write it into Mat_feat_N
        Mat_feat_1(:,i) = reshape(W_n,[],1);
    end
    
    % feature window operation for features of I2
    for i= 1:length(Ftp2(1,:))
        
        % coordinates of current feature
        x = Ftp2(1,i);
        y = Ftp2(2,i);
        
        % window around current feature
        W = I2(y-dist:y+dist,x-dist:x+dist);
        
        % mean, std of window and the normalized Window W_n
        W_mean = (sum(W,"all")/N)*ones(window_length);
        W_std = sqrt((1/(N-1))*trace((W-W_mean)'*(W-W_mean)));
        W_n = (1/W_std)*(W-W_mean);
        
        % convert normalized Window into col-vector and write it into Mat_feat_N
        Mat_feat_2(:,i) = reshape(W_n,[],1);
    end
    
    %% NCC calculations
    % cols of Mat_feat_1 = normalized window vectors of features
    % rows of Mat_feat_2' = normalized window vectors of features
    N = length(Mat_feat_1(:,1));
    NCC_matrix = (1/(N-1))*(Mat_feat_2' * Mat_feat_1); % equals (1/(N-1))*trace(Mf2' * Mf1)
    NCC_matrix(NCC_matrix < min_corr) = 0;  % values smaller min_corr to zero
    
    % find coordinates and values of all non-zero elements
    [y,x,val] = find(NCC_matrix);
    
    % sort rows in descending order of values
    NCC_sorted = sortrows([y,x,val],3,'descend');
    
    % put index of sorted values in a column vector by calculating: y + y_dim * (x-1)
    dim_NCC = size(NCC_matrix);
    sorted_index = [NCC_sorted(:,1) + dim_NCC(1) .* (NCC_sorted(:,2)-1) ];

    %% Correspondeces
    % initialize variables and get dim
    cor = zeros(4,length(sorted_index));
    dim_corr = size(NCC_matrix);
    c = 0;
    
    % for-loop to iterate through correlations in decreasing order
    for i = 1:length(sorted_index)
        
        % get Coordinates of current value in NCC_matrix
        x = ceil(sorted_index(i)/dim_corr(1));      % x-coord in NCC_matrix equals the N-th feature in Ftp1
        y = sorted_index(i) - (x-1) * dim_corr(1);  % y-coord in NCC_matrix equals the N-th feature in Ftp2
        
        % check if feature from Ftp1 was already mapped to a feature in Ftp2
        if NCC_matrix(y,x) ~= 0
            c = c+1;                                % count up number of valid mapped features
            cor(:,c) = [Ftp1(:,x); Ftp2(:,y)];      % write Coordinates of mapped features in cor
            NCC_matrix(:,x) = 0;                    % delete all correlations from curr feature of Ftp1 with features of Ftp2
        end
    end
    cor = cor(:,1:c);                               % cutoff rest of cor with only 0-columns

    %% Visualize the correspoinding image point pairs
    dim = size(cor);
    if do_plot 
        h = imshow(I1);
        hold on
        h.AlphaData = 0.5;
        plot(cor(1,:),cor(2,:),'*red');
        %text(cor(1,:),cor(2,:),string([1:dim(2)]),'Color','green','FontSize',8);
        k = imshow(I2);
        k.AlphaData = 0.5;
        plot(cor(3,:),cor(4,:),'*green');
        for i = 1:length(cor(1,:))
            plot([cor(1,i),cor(3,i)],[cor(2,i),cor(4,i)],'-','color',[0.3010 0.7450 0.9330]);
        end
        axis off
        hold off
    end
end