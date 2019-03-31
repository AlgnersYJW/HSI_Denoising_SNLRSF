
function [Y_rec, nnz_sigular, U_group,singularVals] = LowRankRecovery3(Y, sigma)
 % Low rank matrix recovery
%Input
%Y: Noisy matrix of size B*N
%sigma: noise stand deviation

%Output
%Y_rec: a low rank recovery of Y
%nnz_sigular: Number of nonzero singulars
%U_group: unitary matrices U_group so that Y = U_group*S*V';
%singularVals: filtered singular values

[B, N] = size(Y);

%  [w Rw] = estNoise(Y,noise_type);
% % [bands E]=hysime(Y,w,Rw);    %bands = Ek+2;
% sigma = sqrt(mean(diag(Rw)));

Y = Y./( sqrt(N)* sigma);  

[U,S,V] = svd(Y);

beta = B/N; 
 
lamda = diag(S);
lamda_x = 0.5*( lamda.^2 -1-beta+ sqrt( abs(( lamda.^2 -1 - beta ).^2 - 4*beta ) ) );
lamda_x = sqrt(lamda_x);
 
theta = (1- beta./( lamda_x.^4)) ./ ( 1 + beta./( lamda_x.^2) );
theta = sqrt(theta);
 
phi = (1- beta./( lamda_x.^4)) ./ ( 1 + 1./( lamda_x.^2) );
phi = sqrt(phi);

eta = lamda_x.*theta.*phi;
 
[m,n]=size(S);
eta_dia=zeros(m*n,1); %diagonal matrix
bb=find(lamda  > 1+sqrt(beta));
nnz_sigular = size(bb,1); %Number of nonzero elements
eta_dia((bb-1)*m+bb)=eta(bb);
eta_dia=reshape(eta_dia,[m n]);

Y_rec = U*eta_dia*V' ;

Y_rec = Y_rec.*( sqrt(N)* sigma);
U_group = U ;
singularVals = diag(eta_dia);
end
 