%Inclass 14

%Work with the image stemcells_dapi.tif in this folder

% (1) Make a binary mask by thresholding as best you can

img = imread('stemcells_dapi.tif');
imshow(img, []);

img_mask = img > quantile(quantile(img, 0.8), 0.8);
imshow(img_mask, []);

% (2) Try to separate touching objects using watershed. Use two different
% ways to define the basins. (A) With erosion of the mask (B) with a
% distance transform. Which works better in this case?

CC = bwconncomp(img_mask);
stats = regionprops(CC, 'Area');
area = [stats.Area];
fusedCandidates = area > mean(area) + std(area);
sublist = CC.PixelIdxList(fusedCandidates);
sublist = cat(1, sublist{:});
fusedMask = false(size(img_mask));
fusedMask(sublist) = 1;
imshow(fusedMask, 'InitialMagnification', 'fit');

% eroding
s = round(1.2*sqrt(mean(area))/pi);
nucmin = imerode(fusedMask, strel('disk',s));
imshow(nucmin, 'InitialMagnification', 'fit');

% get outside region
outside = ~imdilate(fusedMask, strel('disk',1));
imshow(outside, 'InitialMagnification', 'fit');

% basins for ws
basin = imcomplement(bwdist(outside));
basin = imimposemin(basin, nucmin | outside);
pcolor(basin);
shading flat;

L = watershed(basin);
imshow(L);
colormap('jet');
caxis([0 20]);

% combining
newNuclearMask = L > 1 | (img_mask - fusedMask);
imshow(newNuclearMask, 'InitialMagnification', 'fit');


% dist transform
dist_t = bwdist(img_mask);
imshow(dist_t,[]);

D = bwdist(~img_mask);
imshow(D, []);

D = -bwdist(~img_mask);
imshow(D, []);

L = watershed(D);
imshow(L, []);

img_mask(L == 0) = 0;
imshow(img_mask);

% the erosion method seems to do a better job in separating the nucleus
% using the distance transform there is an over fragmentation
% separating nuclei of single cells.
% Erosion better identify and separates cells that seem to be merged.
