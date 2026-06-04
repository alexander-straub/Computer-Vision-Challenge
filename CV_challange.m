abs_path = "C:\Users\alex9\Documents\CV_Challange\kicker\images\dslr_images_undistorted\";
%image1 = imread('C:\Users\alex9\Documents\CV_Challange\delivery_area\images\dslr_images_undistorted\DSC_0675.jpg');
%image2 = imread('C:\Users\alex9\Documents\CV_Challange\delivery_area\images\dslr_images_undistorted\DSC_0676.jpg');
numbers = [4,5];
images = load_and_display('\Users\Kilian\MATLAB\Projects\CV_Challange\resources\project\data\example_kicker\images',false,numbers);
img1 = images{4};
img2 = images{5};
gray_image1 = im2gray(img1);
gray_image2 = im2gray(img2);
points1 = detectHarrisFeatures(gray_image1);
strongest = selectStrongest(points1,50); 
points2 = detectHarrisFeatures(gray_image2);
strongest2 = selectStrongest(points2,50);
%correspondences_robust = matchFeatures(points1.Location, points2.Location);
%correspondences_robust = matchFeatures(points1.Location, points2.Location);

[features1,valid_points1] = extractFeatures(gray_image1,points1);
[features2,valid_points2] = extractFeatures(gray_image2,points2);
correspondences_robust = matchFeatures(features1, features2);
matchedpoints1 = valid_points1(correspondences_robust(:,1),:);
matchedpoints2 = valid_points2(correspondences_robust(:,2),:);

figure;
showMatchedFeatures(gray_image1,gray_image2,matchedpoints1,matchedpoints2);

Hs = zeros(3*4,3);
% 4 point algorithm
dim_mp = size(matchedpoints1.Location);
for k = 1:4
    index = ceil(rand(4,1)*dim_mp(1));
    Hs(3*k-2:3*k,:) = four_point_alg(matchedpoints1(index,:), matchedpoints2(index,:), K_matrix);
end


K = cameraIntrinsics([3408.59 3408.87],[3117.24 2064.07],[6208 4135]);
K_matrix = [3408.59 0 3117.24; 0 3408.87 2064.07; 0 0 1.0];

[E,inliers] = estimateEssentialMatrix(matchedpoints1,matchedpoints2,K);
disp(E);
[T1, R1, T2, R2, U, V] = TR_from_E(E);

x1 = matchedpoints1.Location(inliers>0,:);
x2 = matchedpoints2.Location(inliers>0,:);

[T, R, lambda, P1] = reconstruction(T1, T2, R1, R2, [x1'; x2'], K_matrix);
[repro_error, x2_repro] = backprojection([x1'; x2'], P1, gray_image2, T, R, K_matrix);



figure;
h = imshow(gray_image1);
hold on
h.AlphaData = 0.5;
plot(points1.Location(:,1),points1.Location(:,2),'*blue');
k = imshow(gray_image2);
k.AlphaData = 0.5;
plot(points2.Location(:,1),points2.Location(:,2),'*green');
for i = 1:length(correspondences_robust(:,1))
    index = correspondences_robust(i,1);
    index2 = correspondences_robust(i,2);
    plot([points1.Location(index,1),points2.Location(index2,1)],[points1.Location(index,2),points2.Location(index2,2)],'-','color',[0.3010 0.7450 0.9330]);
end
axis off
hold off





