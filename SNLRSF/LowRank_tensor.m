
function  [group_est_weighted, group_weight] = LowRank_tensor( group,idx_similar_cube, sigma)
%Low-rank approximation based filtering for grouped similar patches
%Input
% group: noisy 3D tensor
% idx_similar_cube: index of the location of similar patches
% sigma: noise stand deviation

%Output
% group_est_weighted: denoised 3D tensor
% group_weight = 1, meaning we are using equal weight for estimates of
% a same patch
addpath('tensor_toolbox'); 
if exist('tenmat.m','file') ==0
           errordlg({'Tensor Toolbox not found! ','Download from http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html and install it in the folder .../tensor_toolbox'});
     error('Tensor Toolbox not found! Download from http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html and install it in the folder .../tensor_toolbox');
end

% [ length_pat_square ,N2, b] = size(group);

group_tensor = tensor(group);
for i_mode =   3:-1:1
    
    if i_mode == 3
        group_mode  = tenmat(group_tensor,i_mode);
        group_mode3 = group_mode;
        Y = double(group_mode);
        [ ~, idx, ic] = unique(idx_similar_cube)  ;
        
        [Y_unique_apx, nnz_sigular(i_mode,1), U{i_mode},singularVals{i_mode}] = LowRankRecovery3(Y(:,idx'),sigma);
 
        Y_unique_apx =  U{i_mode}'*Y(:,idx');
        Y_apx = Y_unique_apx(:,ic');
        group_mode(:,:)  = Y_apx(:,:) ;
        group_tensor=tensor(group_mode); %update group_tensor
        group_tensor = double(group_tensor);
        group_tensor(:,:,nnz_sigular(i_mode,1)+1:end) = [];
        group_tensor = tensor(group_tensor);
    else
        group_mode  = tenmat(group_tensor,i_mode);
        Y_input = double(group_mode);
        
        [Y_apx, nnz_sigular(i_mode,1), U{i_mode},singularVals{i_mode}] = LowRankRecovery3(Y_input,sigma);
        if i_mode==1
            Y_apx_stored_mode1 = Y_apx;
        end
       
        Y_apx =  U{i_mode}'*Y_input;
        
        group_mode(:,:)  = Y_apx(:,:) ;
        group_tensor=tensor(group_mode);
        group_tensor = double(group_tensor);
        switch i_mode
            case 2
                group_tensor(:, nnz_sigular(i_mode,1)+1:end,:) = [];
                group_mode2 = group_mode;
            case 1
                group_tensor( nnz_sigular(i_mode,1)+1:end,:,:) = [];
                group_mode1 = group_mode;
                
        end
        group_tensor = tensor(group_tensor);
    end
    
end


i_mode = 1;
if 0
    group_mode = tenmat(group_tensor,i_mode);
    Y_mode = double(group_mode);
    Y_mode = U{i_mode}(:,1:nnz_sigular(i_mode,1))*Y_mode;
    group_mode1(:,:) = Y_mode(:,:);
else
    group_mode1(:,:) = Y_apx_stored_mode1(:,:);
end
group_mode1_tensor =  tensor( group_mode1);

i_mode =   2;
Y_mode = tenmat(group_mode1_tensor,i_mode);
Y_mode = double(Y_mode);
Y_mode = U{i_mode}(:,1:nnz_sigular(i_mode,1))*Y_mode;
group_mode2(:,:) = Y_mode(:,:);
group_mode2_tensor =  tensor( group_mode2);

i_mode =   3;
if nnz_sigular(i_mode,1)>1
    Y_mode = tenmat(group_mode2_tensor,i_mode);
    Y_mode = double(Y_mode);
else
    tmp = double(group_mode2_tensor);
    Y_mode = tmp(:)';
end
Y_mode = U{i_mode}(:,1:nnz_sigular(i_mode,1))*Y_mode;
group_mode3(:,:) = Y_mode(:,:);
group_mode3_tensor =  tensor( group_mode3);

group_est_weighted =  double(group_mode3_tensor);

group_weight = 1 ;

end