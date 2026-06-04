%% Load images

    %images = load_and_display('\\10.1.1.11\kilian\Documents\Studium\Matlab\resources\project\Code\data\delivery_area_dslr_undiplottingstorted\delivery_area\images\dslr_images_undistorted',false);
    images = load_and_display('\Users\Kilian\MATLAB\Projects\CV_Challange\resources\project\data\delivery_area_dslr_undistorted\delivery_area\images\dslr_images_undistorted',false,[1,2]);
    %images = load_and_display('\Users\Kilian\MATLAB\Projects\CV_Challange\resources\project\data\example_kicker\images',false,numbers);
    img1 = images{1};
    img2 = images{2};
    gray1 = rgb_to_gray(img1);
    gray2 = rgb_to_gray(img2);
    
    imshow(gray1)
    imshow(gray2)
    % test

%% Calculate Harris features
    min_dist = 20;
    N = 40;
    tile_size = [200, 300];
    [features1,~,~] = harris_detector(gray1, 'min_dist', min_dist, 'N', N, 'tile_size', tile_size, 'do_plot', false);
    [features2,~,~] = harris_detector(gray2, 'min_dist', min_dist, 'N', N, 'tile_size', tile_size, 'do_plot', false);

%% Correspondence estimation
    correspondences = point_correspondence(gray1,gray2,features1,features2,'window_length',55,'min_corr', 0.9);

%% Determine robust corresponding image points with the RANSAC algorithm
    correspondences_robust = F_ransac(correspondences, 'tolerance', 0.04);   %% 0.04

%% Visualize robust corresponding image points
    figure;
    h = imshow(gray1);
    hold on
    h.AlphaData = 0.5;
    plot(correspondences(1,:),correspondences(2,:),'*red');
    k = imshow(gray2);
    k.AlphaData = 0.5;
    plot(correspondences(3,:),correspondences(4,:),'*green');
    for i = 1:length(correspondences(1,:))
    plot([correspondences(1,i),correspondences(3,i)],[correspondences(2,i),correspondences(4,i)],'-','color',[0.3010 0.7450 0.9330]);
    end
    axis off
    hold off

    figure;
    h = imshow(gray1);
    hold on
    h.AlphaData = 0.5;
    plot(correspondences_robust(1,:),correspondences_robust(2,:),'*red');
    k = imshow(gray2);
    k.AlphaData = 0.5;
    plot(correspondences_robust(3,:),correspondences_robust(4,:),'*green');
    for i = 1:length(correspondences_robust(1,:))
    plot([correspondences_robust(1,i),correspondences_robust(3,i)],[correspondences_robust(2,i),correspondences_robust(4,i)],'-','color',[0.3010 0.7450 0.9330]);
    end
    axis off
    hold off

%% Calculate essential matrix
K = [3408.59 0 3117.24; 0 3408.87 2064.07; 0 0 1];
E = epa(correspondences_robust, K);
disp(E);

%% Draw epipolarlines
dim = size(correspondences_robust);
homo_cal_corr_rob = [inv(K)*[correspondences_robust(1:2,:); ones(1,dim(2))]; inv(K)*[correspondences_robust(3:4,:); ones(1,dim(2))]];
Rot = [cos(pi), sin(pi), 0; -sin(pi), cos(pi), 0; 0, 0, 1];
l1 = E'*homo_cal_corr_rob(4:6,:);
l2 = E*homo_cal_corr_rob(1:3,:);
%l1 = E'*correspondences_robust(4:6,:);
%l2 = E*correspondences_robust(1:3,:);


lines = [correspondences_robust(1:2,:); correspondences_robust(1:2,:) + 10000*[l1(1,:);-l1(2,:)]; correspondences_robust(3:4,:); correspondences_robust(3:4,:) + 10000*[l2(1,:);-l2(2,:)]];

figure;
subplot(2,1,1);
imshow(gray1);
hold on
for i = 1:length(lines(1,:))
    plot([lines(1,i),lines(3,i)],[lines(2,i),lines(4,i)],'-','color',[0.3010 0.7450 0.9330]);
end 
subplot(2,1,2);
imshow(gray2);
hold on
for i = 1:length(lines(1,:))
    plot([lines(5,i),lines(7,i)],[lines(6,i),lines(8,i)],'-','color',[0.3010 0.7450 0.9330]);
end 

%% Determine R and T
[T1, R1, T2, R2, U, V] = TR_from_E(E);
[T, R, lambda, P1] = TR_calculation(T1, T2, R1, R2, correspondences_robust, K);
[repro_error, x2_repro] = backprojection(correspondences_robust, P1, img2, T, R, K);

%% Find points in planars
k = 20;
indexes = zeros(4,k*20);
vs = zeros(3*k,4);
Hs = zeros(3*k,3);
for i = 1:k
    index = ceil(rand(4,1)*dim(2));
    indexes(:,i) = index;
    Hs(3*i-2:3*i,:) = four_point_alg_ess(homo_cal_corr_rob(:,index), K);
    for j = 1:4
        vs(3*i-2:3*i,j) = three_point_system(E,T,homo_cal_corr_rob(1:3,index(j)),homo_cal_corr_rob(4:6,index(j)));

    end
end
