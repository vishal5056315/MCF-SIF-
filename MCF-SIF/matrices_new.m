function metrics = matrices_new(xrcw,a,b)
% MATRICES_NEW Compute fusion metrics and Petrovic-based QABF/LABF/NABF
% Inputs:
%   xrcw - fused image (uint8 or double)
%   x    - cell array with source images: x{1}, x{2}
% Output:
%   metrics - struct containing computed metrics (API, SD, AG, Entropy, MIF, FS1, corr, SF, QABF, LABF, NABF, NABF1)

% Ensure images are double [0,1] for computations where needed
x1 = im2double(a);
x2 = im2double(b);
% For histogram/joint-hist operations, use uint8 0..255
F_uint8 = img_to_uint8(xrcw);
A_uint8 = img_to_uint8(x1);
B_uint8 = img_to_uint8(x2);

% Use double version for numerical computations
F = im2double(F_uint8);
[p,q] = size(F);

%% API
API=mean(F_uint8(:))

%%% Standard Deviation (SD).
[p,q]=size(F_uint8);
xrcwd=double(F_uint8);
SD=sqrt(sum(sum((xrcwd-API).^2))/(p*q))

%%% Average Gradient (AG).
tmp1=[xrcwd(2:p,:);xrcwd(1,:)];
tmp2=[xrcwd(:,2:q) xrcwd(:,1)];
tmp=sqrt((xrcwd-tmp1).^2+(xrcwd-tmp2).^2);
AG=sum(tmp(:))/(p*q)
clear tmp1 tmp2 tmp

%% Entropy
histF = imhist_fn(F_uint8);
pdfF = histF / (p*q);
entropyF = -sum(pdfF(pdfF>0).*log2(pdfF(pdfF>0)));
metrics.Entropy = entropyF;
fprintf('Entropy = %.6f\n', entropyF);

%% Mutual Information (MIF) between each source and fused image
histA = imhist_fn(A_uint8); pdfA = histA / (p*q);
histB = imhist_fn(B_uint8); pdfB = histB / (p*q);

jpdfAF = joint_hist_fn(A_uint8, F_uint8) / (p*q);
jpdfBF = joint_hist_fn(B_uint8, F_uint8) / (p*q);
jpdfAB = joint_hist_fn(A_uint8, B_uint8) / (p*q);

% MIAF and MIBF: compute mutual info using matrix operations
MIAF = mutual_info_from_joint(jpdfAF, pdfA, pdfF);
MIBF = mutual_info_from_joint(jpdfBF, pdfB, pdfF);
MIF = MIAF + MIBF;
metrics.MIAF = MIAF;
metrics.MIBF = MIBF;
metrics.MIF = MIF;
fprintf('MIAF = %.6f, MIBF = %.6f, MIF = %.6f\n', MIAF, MIBF, MIF);

%% Fusion Symmetry (FS1) (as in your code)
if (MIAF + MIBF) == 0
    FS1 = NaN;
else
    FS1 = 2 - abs((MIAF/(MIAF+MIBF)) - 0.5);
end
metrics.FS1 = FS1;
fprintf('FS1 = %.6f\n', FS1);

%% Average Normalized Correlation (rAF, rBF, corr)
diffF = F - mean(F(:));
diffA = x1 - mean(x1(:));
diffB = x2 - mean(x2(:));
rAF = sum(diffA(:).*diffF(:)) / sqrt(sum(diffA(:).^2)*sum(diffF(:).^2) + eps);
rBF = sum(diffB(:).*diffF(:)) / sqrt(sum(diffB(:).^2)*sum(diffF(:).^2) + eps);
corr_avg = (rAF + rBF)/2;
metrics.rAF = rAF;
metrics.rBF = rBF;
metrics.corr = corr_avg;
fprintf('rAF = %.6f, rBF = %.6f, corr = %.6f\n', rAF, rBF, corr_avg);

%% Spatial Frequency (SF)
xrcwd = double(F_uint8);
rf = sqrt(sum(sum((xrcwd(:,2:q) - xrcwd(:,1:q-1)).^2)) / (p*q));
cf = sqrt(sum(sum((xrcwd(2:p,:) - xrcwd(1:p-1,:)).^2)) / (p*q));
SF = sqrt(rf^2 + cf^2);
metrics.SF = SF;
fprintf('SF = %.6f\n', SF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Petrovic-based 3 Metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters (as in your code)
Td = 2; wt_min = 0.001; P = 1; Lg = 1.5;
Nrg = 0.9999; kg = 19; sigmag = 0.5;
Nra = 0.9995; ka = 22; sigmaa = 0.5;

% Edge Strength & Orientation using sobel
[gvA, ghA] = sobel_fn(A_uint8);
[gvB, ghB] = sobel_fn(B_uint8);
[gvF, ghF] = sobel_fn(F_uint8);
gA = sqrt(gvA.^2 + ghA.^2);
gB = sqrt(gvB.^2 + ghB.^2);
gF = sqrt(gvF.^2 + ghF.^2);

% Preallocate
gAF = zeros(p,q); gBF = zeros(p,q);
aA = zeros(p,q); aB = zeros(p,q); aF = zeros(p,q);

% Compute edge strength ratios and orientations (vectorized where possible)
% Avoid divisions by zero by adding eps where needed.
% Orientation (atan2) is numerically safer than atan(gv/gh)
aA = atan2(gvA, ghA);
aB = atan2(gvB, ghB);
aF = atan2(gvF, ghF);

% gAF and gBF (ratio rule from your original loops)
% Use elementwise logic:
maskA = (gA==0) | (gF==0);
gAF(~maskA) = min(gA(~maskA)./gF(~maskA), gF(~maskA)./gA(~maskA));
maskB = (gB==0) | (gF==0);
gBF(~maskB) = min(gB(~maskB)./gF(~maskB), gF(~maskB)./gB(~maskB));

% aAF, aBF (orientation difference transform)
aAF = abs(abs(aA - aF) - pi/2) * 2/pi;
aBF = abs(abs(aB - aF) - pi/2) * 2/pi;

% Edge Preservation Coefficient
QgAF = Nrg ./ (1 + exp(-kg*(gAF - sigmag)));
QaAF = Nra ./ (1 + exp(-ka*(aAF - sigmaa)));
QAF = sqrt(QgAF .* QaAF);

QgBF = Nrg ./ (1 + exp(-kg*(gBF - sigmag)));
QaBF = Nra ./ (1 + exp(-ka*(aBF - sigmaa)));
QBF = sqrt(QgBF .* QaBF);

% Weight maps wtA and wtB
wtA = wt_min * ones(p,q);
wtB = wt_min * ones(p,q);
mask_gt = (gA >= Td);
wtA(mask_gt) = gA(mask_gt).^Lg;
mask_gtB = (gB >= Td);
wtB(mask_gtB) = gB(mask_gtB).^Lg;

wt_sum = sum(wtA(:) + wtB(:));

% QABF
QABF = (sum(QAF(:).*wtA(:)) + sum(QBF(:).*wtB(:))) / (wt_sum + eps);
metrics.QABF = QABF;
fprintf('QABF = %.6f\n', QABF);

% LABF (fusion loss)
rr = (gF <= gA) | (gF <= gB);
LABF = sum(sum(rr .* ((1 - QAF).*wtA + (1 - QBF).*wtB))) / (wt_sum + eps);
metrics.LABF = LABF;
fprintf('LABF = %.6f\n', LABF);

% NABF1 (Petrovic)
na1 = zeros(p,q);
mask_na1 = (gF > gA) & (gF > gB);
na1(mask_na1) = 2 - QAF(mask_na1) - QBF(mask_na1);
NABF1 = sum(na1(:) .* (wtA(:) + wtB(:))) / (wt_sum + eps);
metrics.NABF1 = NABF1;
fprintf('NABF1 = %.6f\n', NABF1);

% NABF (Kumar)
na = zeros(p,q);
na(mask_na1) = 1;
NABF = sum(sum(na .* ((1 - QAF).*wtA + (1 - QBF).*wtB))) / (wt_sum + eps);
metrics.NABF = NABF;
fprintf('NABF = %.6f\n', NABF);

end

%% ------------------------ Helper local functions ------------------------
function out = img_to_uint8(img)
% Convert image to uint8 [0..255] robustly
if isfloat(img)
    % If already in 0..1, scale; if outside, clip then scale
    minv = min(img(:)); maxv = max(img(:));
    if maxv <= 1 && minv >= 0
        out = im2uint8(img);
    else
        % normalize to 0..1 then uint8
        out = im2uint8(mat2gray(img));
    end
elseif isinteger(img)
    % convert to uint8 preserving scale if possible
    out = im2uint8(img);
else
    out = uint8(img);
end
end

function h = imhist_fn(img_uint8)
% Compute histogram (256 bins) for uint8 image
if ~isa(img_uint8,'uint8')
    img_uint8 = img_to_uint8(img_uint8);
end
h = histcounts(img_uint8(:), 0:256)'; % 256x1
end

function J = joint_hist_fn(I_uint8, J_uint8)
% Joint histogram (256x256) for two uint8 images of same size
if ~isa(I_uint8,'uint8'), I_uint8 = img_to_uint8(I_uint8); end
if ~isa(J_uint8,'uint8'), J_uint8 = img_to_uint8(J_uint8); end
assert(isequal(size(I_uint8), size(J_uint8)), 'Inputs must be same size for joint histogram.');

% Linear indexing to build joint histogram efficiently
Iv = double(I_uint8(:));   % 0..255
Jv = double(J_uint8(:));   % 0..255
idx = Iv*256 + Jv + 1;     % map pair (a,b) to unique index (1..65536)
Jlin = accumarray(idx, 1, [256*256, 1]);
J = reshape(Jlin, [256,256]);
end

function MI = mutual_info_from_joint(jpdf, pdfX, pdfY)
    [rows, cols] = size(jpdf);
    MI = 0;
    for i = 1:rows
        for j = 1:cols
            if jpdf(i,j) > 0 && pdfX(i) > 0 && pdfY(j) > 0
                MI = MI + jpdf(i,j) * log2(jpdf(i,j) / (pdfX(i) * pdfY(j)));
            end
        end
    end
end


function [gv, gh] = sobel_fn(img_uint8)
% Compute vertical and horizontal sobel responses (double)
if ~isa(img_uint8,'uint8'), img_uint8 = img_to_uint8(img_uint8); end
I = double(img_uint8);
hx = [1 0 -1; 2 0 -2; 1 0 -1]; % horizontal kernel (gh)
hy = hx';                    % vertical kernel (gv)
gh = conv2(I, hx, 'same');   % horizontal gradient
gv = conv2(I, hy, 'same');   % vertical gradient
end
