clear all;
close all;
clc;

% Load image
I = imread("lake_fused.png");

% Convert to grayscale if RGB
if size(I, 3) == 3
    I = rgb2gray(I);
end

% Extract non-zero pixels only
nonzero_pixels = double(I(I > 0));

% Compute statistics
mean_val   = mean(nonzero_pixels);
median_val = median(nonzero_pixels);
mode_val   = mode(nonzero_pixels);
std_val    = std2(I);   % now should work
var_val    = var(nonzero_pixels);

% Display results
disp(['Mean: ', num2str(mean_val)]);
disp(['Median: ', num2str(median_val)]);
disp(['Mode: ', num2str(mode_val)]);
disp(['Standard Deviation: ', num2str(std_val)]);
disp(['Variance: ', num2str(var_val)]);