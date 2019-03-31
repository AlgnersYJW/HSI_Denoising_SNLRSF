%==========================================================================
% Using the code should cite the following paper
% 
% C. Cao, J. Yu, C. Zhou, K. Hu, F. Xiao and X. Gao, "Hyperspectral Image Denoising via Subspace-Based 
% Nonlocal Low-Rank and Sparse Factorization," in IEEE Journal of Selected Topics in 
% Applied Earth Observations and Remote Sensing, vol. 12, no. 3, pp. 973-988, March 2019.
%
% by Jie Yu, 2018.
%==========================================================================

%% HSI Denoising Simulated Experiment 
clc; clear all; close all;
addpath(genpath(pwd));

%% Load simulated noise HSI data
% Please first running Demo_Generate_Simulated_Data
load Noisy_WDC_CASE1.mat;
% load Noisy_PaviaU_Case1.mat;
[nr,nc,L] = size(Noisy_Img);

X = Img;
Y = Noisy_Img;
clear Img Noisy_Img;

%% Subspace NonLocal Low Rank and Sparse Factorization
t1=clock;
sub_dim = 5;    % selected by different dataset
[ SNLRSF_Ys ] = SNLRSF_HSI_Denoising( Y, sub_dim );
t2=clock;
time = etime(t2,t1);
[mpsnr,psnr] = MPSNR(X,SNLRSF_Ys);
[mssim,ssim] = MSSIM(X,SNLRSF_Ys);
ergas        = ErrRelGlobAdimSyn(X,SNLRSF_Ys);
msa          = MSA(X, SNLRSF_Ys);
% save WDC_Case1_SNLRSF SNLRSF_Ys mpsnr psnr mssim ssim ergas msa time