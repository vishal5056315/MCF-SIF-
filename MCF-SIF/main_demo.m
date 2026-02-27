% MCF-SIF: Mean Curvature Filter-based Structure-aware Image Fusion
%
% This script is the main file for running the experiments from:
%   <Author Names:Vishal, Vinay Kukreja, Ayush Dogra, Bhawna Goyal>
%   "Photoacoustic Biomedical Signal Fusion for Multi-Patch Based Vascular Details: A Curvature-Based Approach"
%
% It demonstrates:
%   - Loading multi-patch photoacoustic vascular images (rat ear).
%   - Running the MCF-SIF fusion algorithm implemented in MCF_SIF.m.
%   - Computing quantitative metrics using matrices_new.m.
%
% Structure-aware fusion components (saliency/decision maps and
% structure-preserving behavior) are related to:
%   <Author Names: Wen Li, Yuange Xie, Haole Zhou, Ying Han, and Kun Zhan>
%   "Structure-aware image fusion" 
%   <Journal: Optik> 
%   <Year:2018>.
%
% mean_curvature_filter.m
% -------------------------------------------------------------------------
% Iterative mean curvature motion PDE for geometry-aware weight refinement.
%
% This implementation follows the standard mean curvature motion formulation
% used for anisotropic smoothing and edge-preserving evolution of level sets.
%
% If you reuse this function in your own work, please cite:
%   <Author Names: Yuanhao Gong and Ivo F. Sbalzarini>,
%   "<Curvature filters efficiently reduce certain variational energies>,"
%   <Journal: IEEE Transactions on Image Processing>
%   <Year:2017>.
%
% -------------------------------------------------------------------------

% If you use this script or the associated fusion method in your work,
% please cite both our paper and the relevant baseline references.
%
% -------------------------------------------------------------------------

clear;
close all;clc
addpath(genpath(pwd));

I1 = imread('lake\1.png');
I2 = imread('lake\2.png');

figure; imshow(I1); title("Input 1");
figure; imshow(I2); title("Input 2");

% ---------------------------------------------------------
% Ensure both images are RGB (convert if grayscale)
% ---------------------------------------------------------
if size(I1,3) == 1
    I1 = repmat(I1, [1 1 3]);
end
if size(I2,3) == 1
    I2 = repmat(I2, [1 1 3]);
end
tic
F_rgb = zeros(size(I1), 'uint8');

for c = 1:3
    F_rgb(:,:,c) = MCF_SIF(I1(:,:,c), I2(:,:,c));
end

toc

figure; imshow(F_rgb); title("Fused Image");
imshow(F_rgb)
imwrite(F_rgb,"S1.png");
matrices_new(rgb2gray(F_rgb),im2gray(I1),im2gray(I2));
