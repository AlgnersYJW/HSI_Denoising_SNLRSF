
function D = SNLRSF_HSI_Denoising(Y, sub_dim)

[row,col,band] = size(Y);

YD = reshape(Y, row*col, band)';
noise_type='additive';  verbose1 = 'off';
[w, Rw]    = estNoise_snlrsf(YD,noise_type,verbose1);
Rw_save = Rw;
YD = sqrt(inv(Rw_save))*YD;

% Estimate the orthogonal matrix E by HySime;
[w, Rw]    = estNoise_snlrsf(YD,noise_type,verbose1);
[kf,E_all] = hysime_snlrsf(YD,w,Rw,verbose1);
if ~strcmp(verbose1,'on'),fprintf(1,'The signal subspace dimension is: k = %d\n',kf);end
E = E_all(:,1:sub_dim); % size of band*sub_dim

% parameters for nonlocal low rank recovery
t          = 39;  % radius of window for searching similar patches, 
                  % meaning size of the search windiow is (2*t+1)*(2*t+1)
f          = 5;   % length of a patch, meaning size of a patch is f*f.
NL_num     = 8;  % number of similar blocks
refer_step = 3;   % sliding step to process every next reference patch

% regularized paramaters
lam    = 0.1;   % for low rank term  lambda1
tau    = 0.05;  % for sparse term  lambda2

AL_iters = 100;
iter = 1;

% initilizations
Y = YD; clear YD

Zt = E'*Y;
Z = reshape(Zt', row, col, sub_dim);
% update L
[Z, Spa_Img, Spa_W] = NonLocal_LR(Z, t, f, NL_num, refer_step);
%     load L
Z = reshape(Z, row*col, sub_dim)';
Spa_Wei = repmat(Spa_W, sub_dim, 1);
while iter <= AL_iters    
    % update S
    S = MySoftTh( Y-E*Z, tau );    
    % update Z
    Z = (lam*Spa_Img + E'*(Y-S)) ./ (lam*Spa_Wei + 1);    
    % update E
    E_est   = (Y-S)*Z';
    [U,~,V] = svd(E_est,'econ');
    E = U*V';    
    iter = iter + 1;
end

D = E*Z;
D = sqrt(Rw_save)*D;
D = reshape(D', row, col, band);
end

%% This is soft thresholding operation
function X= MySoftTh(B,lambda)
X=sign(B).*max(0,abs(B)-(lambda));
end

%% nonlocal low rank recovery
function [output_img, Y_denoised, count_Y_denoised] = NonLocal_LR(img_noisy, t, f, N2, step)

[r, c, b] = size(img_noisy);
N=r*c;
Y_noisy = reshape(img_noisy, N, b)';

% Extract features : SVD
% apply SVD on img_noisy, then use first p components to compute similarity between two patches.
p   = 1;  
[U_ss,D]=svd(Y_noisy,'econ');
U_ss(:,p+1:end) = []; % remain first p component
Y_feature =U_ss'*Y_noisy; % truncated SVD
img_feature = reshape(Y_feature', r, c, p);

length_pat = f;

% Replicate the boundaries of the input image
img_feature2 = padarray(img_feature,[f-1 f-1],'symmetric','post');

% total patches for the image
for ip = 1:p
    T_pat(:,:,ip) = im2col( img_feature2(:,:,ip), [length_pat,length_pat], 'sliding');
end

%index image
idx_img = reshape( 1: N , [r, c]);
idx_img = padarray( idx_img, [f-1 f-1],'symmetric','post');
idx_img_pat=im2col( idx_img, [length_pat,length_pat], 'sliding');

kernel = ones(length_pat^2,1)/length_pat^2; % each entry in a patch is equally important.

i_patch_count=0;
Y_denoised = zeros(size( Y_noisy));
count_Y_denoised = zeros(1,size(Y_noisy,2));
for j = 1:step:c % col 
    for i = 1:step:r % row    
        i_patch_count = i_patch_count+1;
        i1 = i;  j1 = j;
        
        % reference patch id
        ref_pat_idx = idx_img(i1, j1);
        ref_pat = T_pat(:,ref_pat_idx,:);
        
        % define four cordinate of search window
        rmin = max(i1-t,1);
        rmax = min(i1+t,r+f-1);
        smin = max(j1-t,1);
        smax = min(j1+t,c+f-1);
        
        % other patchs in the search window
        patches_idx = idx_img(rmin:1:rmax,smin:1:smax);
        patches_idx = patches_idx(:);
        patches     = T_pat(:,patches_idx,:);
        
        % calculate distance between reference patch and other patchs
        dis =zeros(p,size(patches,2));
        for ip=1:p
            dis(ip,:) = sum( kernel * ones(1,size(patches, 2)).* bsxfun(@minus,patches(:,:,ip),ref_pat(:,1,ip)).^2);  
        end
        
        % sort distance and extract similar patchs index with define number
        w = -sum(dis,1);
        [~, I] = sort(w, 'descend');
        I = I(1: N2);

        if length(I)>= 1
            % non-local similar patchs
            patches_sel_idx     = patches_idx(I,1)';
            idx_similar_cube_2d = idx_img_pat(:,patches_sel_idx);  
            idx_similar_cube    = idx_similar_cube_2d(:)';
            
            % extract similar patch and form group
            group = Y_noisy(:,idx_similar_cube);
            group = reshape(group', length_pat^2,N2, b);
            
            % denoise in each group
            sigma = 1; %noise std of each entry of img_noisy is 1, since noise has been whiten.
            [group_est, weight] = LowRank_tensor( group,idx_similar_cube, sigma);
         
            for ic = 1:N2  % patch by patch
                tmp =  idx_similar_cube_2d(:,ic)';
                Y_denoised(:,tmp) = Y_denoised(:,tmp)+squeeze(group_est(:,ic,:))';
                count_Y_denoised(1,tmp) = count_Y_denoised(1,tmp)+weight;
            end
        end
    end
end

output_img = bsxfun(@rdivide,Y_denoised,count_Y_denoised);
output_img=reshape(output_img',r, c, b);
end
