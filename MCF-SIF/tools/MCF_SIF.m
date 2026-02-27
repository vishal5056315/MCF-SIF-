function F = MCF_SIF(Sa,Sb)
% if ~exist('sigma_r','var')
%     sigma_r = 0.2;
% end

%% normalization
GA = im2double(Sa);
GB = im2double(Sb);
%% Smoothing
r = 3;    lambda = 0.01;
[hei, wid] = size(GA);
N = boxfilter(ones(hei, wid), r);
Ga = smoothing(GA, r, lambda,N);
Gb = smoothing(GB, r, lambda,N);
% r = 7;
N = boxfilter(ones(hei, wid), r);
%% Structure
h = [1 -1];

mc_iters = 2;     % iterations (tune for smoothness)
mc_dt    = 0.01; 
MA = abs(conv2(Ga,h,'same')) + ...
     abs(conv2(Ga,h','same'));
MB = abs(conv2(Gb,h,'same')) + ...
     abs(conv2(Gb,h','same'));
DA = MA - MB;
IA = boxfilter(DA,r) ./ N>0;
        for t = 1:3
            IA = double(IA > 0.5);
            IA = mean_curvature_filter(IA,mc_iters,mc_dt);
        end
%         imwrite(IA,'ia1.bmp');
DB = MB - MA;
IB = boxfilter(DB,r) ./ N>0;
        for t = 1:3
            IB = double(IB > 0.5);
            IB = mean_curvature_filter(IB,mc_iters,mc_dt);
        end
%% Result
F = IA.*GA + (IB).*GB;
F = uint8(255*F);