
%% Demo on generate simulate HSI data
clc; clear all;
addpath('Data');

%% Load HSI
load('Ori_WDC.mat');          % not normalize
% load('Ori_PaviaU.mat');     % not normalize

%% nomalize the image
[nr,nc,L] = size(Img);
%nomorlize the observed X
Img = Img./repmat(max(max(Img,[],1),[],2),nr,nc);

%% Add noise in different ways
Noisy_Img = Img;
CaseType = 1;
switch CaseType
    case 1
        % zero-mean Gaussian noise with the same standard deviation
        sigma = 0.1;
        Noisy_Img = Img + sigma*randn(nr,nc,L);
        save Noisy_WDC_CASE1 Img Noisy_Img sigma;
    case 2
        % zero-mean Gaussian noise with the different standard deviation
        % randomly sampled from [0.1,0.2] was added to the different band
        sig =  0.1+(0.2-0.1)*rand(L,1);
        for i = 1 : L
            Noisy_Img(:,:,i) = Img(:,:,i) + sig(i)*randn(nr,nc);
        end
        sigma = mean(sig);
        save Noisy_WDC_CASE2 Img Noisy_Img sigma;
    case 3
        % Gaussian noise the same as case 2, and impulse noise with a
        % density percentage of 20% was added from band 90 to band 130
        sig =  0.1+(0.2-0.1)*rand(L,1);
        for i = 1 : L
            Noisy_Img(:,:,i) = Img(:,:,i) + sig(i)*randn(nr,nc);
        end
        sigma = mean(sig);
        ratio = 0.2;
        RL_sp = randi(L,20,1);
        for i = 1 : 20
            Noisy_Img(:,:,RL_sp(i)) = imnoise(Noisy_Img(:,:,RL_sp(i)), 'salt & pepper', ratio);
        end
        save Noisy_WDC_CASE3 Img Noisy_Img sigma ratio;
    case 4
        %  Gaussian noise and impulse noise the same as case 3
        sig =  0.1+(0.2-0.1)*rand(L,1);
        for i = 1 : L
            Noisy_Img(:,:,i) = Img(:,:,i) + sig(i)*randn(nr,nc);
        end
        sigma = mean(sig);
        ratio = 0.2;
        RL_sp = randi(L,20,1);
%         RL_sp = sort(RL_sp);
        for i = 1:20
            Noisy_Img(:,:,RL_sp(i)) = imnoise(Noisy_Img(:,:,RL_sp(i)), 'salt & pepper', ratio);
        end
        % dead lines were added from band 126 to band 145 with the number
        % of deadlines being randomly selected from 3 to 10, and the width
        % of each dead line was randomly generated from 1 to 3
        RL_dl = randi(L,10,1);
        RL_dl = [RL_sp(1:10); RL_dl];
%         RL_dl = sort(RL_dl);
        for i=1:20
            indp=randperm(8,1)+2;
            ind=randperm(nc-1,indp);
            an = zeros(length(ind),1);
            for k=1:length(ind)
                an(k)=randperm(3,1);
            end
            % searching the location of an which value is 1,2,3
            loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
            Noisy_Img(:,ind(loc1),RL_dl(i))=1; 
            Noisy_Img(:,ind(loc2):ind(loc2)+1,RL_dl(i))=1;
            Noisy_Img(:,ind(loc3)-1:ind(loc3)+1,RL_dl(i))=1;
        end
        % dead lines by another way
%         for i=126:145
%             loc_dl = ceil(rand(1,20)*nc);
%             loc_dl = [loc_dl, 20:22];
%             loc_dl = [loc_dl, 142:143];
%             Noisy_Img(:,loc_dl,i) = ones(nr,size(loc_dl,2));
%         end  
        save Noisy_WDC_CASE4 Img Noisy_Img RL_sp RL_dl sigma;
end
