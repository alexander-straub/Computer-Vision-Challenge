function imageColl = load_and_display(dataFolder, doplot, numbers)
%dataFolder = '/home/kilian/MATLAB/Projects/CV_Challenge/resources/project/data/kicker/images/dslr_images_undistorted';
if ~isdir(dataFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', dataFolder);
  uiwait(warndlg(errorMessage));
  return;
end
filePattern = fullfile(dataFolder, '*.JPG');
jpegFiles = dir(filePattern);
imageColl = cell(1,length(jpegFiles));
for k = numbers
  baseFileName = jpegFiles(k).name;
  fullFileName = fullfile(dataFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  imageArray = imread(fullFileName);
  imageColl{k} = imageArray;
  % imshow(imageArray);  % Display image.
  % drawnow; % Force display to update immediately.
end
if doplot 
    montage(imageColl);
end

% path for delivery area:   \\robbi\kilian\Documents\Studium\Matlab\resources\project\Code\data\delivery_area_dslr_undistorted\delivery_area\images\dslr_images_undistorted
% path for kicker:          \\robbi\kilian\Documents\Studium\Matlab\resources\project\Code\data\kicker_dslr_undistorted\kicker\images\dslr_images_undistorted
% pc oben                   \\10.1.1.11\kilian\Documents\Studium\Matlab\resources\project\Code\data\kicker_dslr_undistorted\kicker\images\dslr_images_undistorted
end
