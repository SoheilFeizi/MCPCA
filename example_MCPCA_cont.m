% Example script
% MCPCA on continuous data

% Citation:
% Soheil Feizi, David Tse, Maximally Correlated Principal Component Analysis, arXiv:1702.05471

clc
close all
clear all

%*************************************
% Parameters of generaing data
n=1000; % number of samples
p=50; % number of variables
qq=10; % dimension of underlying data

% piece-wise linear transformation
%control_fun=1;
B=10; % transforming data with B degree picewise functions

% polynomial transformation
control_fun=2;

std_fun=100;
noise_scale=1;

%*************************************
%MCPCA parameters

q=qq;
d=10; % functions degree freedom
num_init=10; % number of initializations
num_iter=10; % number of repeats of coordinate block descent updates

%*************************************
% Generate latent data samples
X00=randn(n,qq)*randn(qq,p);
[U00,S00,V00]=svd(X00,0);
[sigma_vec00,I00]=sort(diag(S00),'descend');
V00=V00(:,I00);
U00=U00(:,I00);

% low dimensional embedding
U00_scaled=U00;
for i=1:size(U00,2)
    U00_scaled(:,i)=U00(:,i)*sigma_vec00(i);
end

% adding noise 
noise_mat=randn(n,p).*noise_scale;
X0=X00+noise_mat;
X0=normalize_matrix(X0);

% compute Ky Fan norm
[V0,D0]=eig(X0'*X0/n);
[lambda_vec0,I0]=sort(diag(D0),'descend');
kyfan0=kyfan(lambda_vec0');

%*************************************
% Generate observed data samples
% transforming data

X_t=X0*0;
if control_fun==1
% piece-wise linear transformation
    for i=1:p
        [~,cent_vec]=discrt_v2(X0(:,i),B);
        
        y_cent_vec=zeros(B,1);
        y_cent_vec(1)=exprnd(std_fun);
        for j=2:B
            y_cent_vec(j)=y_cent_vec(j-1)+exprnd(std_fun);
        end
        
        if rand(1)>0.5
            y_cent_vec=y_cent_vec;
        else
            y_cent_vec=-y_cent_vec;
        end
        X_t(:,i)=piecewise_lin(X0(:,i),cent_vec,y_cent_vec);
    end
elseif control_fun==2
% polynomial transformation    
    for i=1:p
        
        if rand(1)<1/3
            X_t(:,i)=X0(:,i);
        elseif rand(1)<2/3
            X_t(:,i)=X0(:,i).^3;
        else
            X_t(:,i)=X0(:,i).^5;
        end       
    end   
end

X_t=normalize_matrix(X_t);

[U_t,S_t,V_t]=svd(X_t);
[sigma_vec_t,I_t]=sort(diag(S_t),'descend');
V_t=V_t(:,I_t);
U_t=U_t(:,I_t);

% low dimensional embedding using PCA
U_t_scaled=U_t;
for i=1:size(U_t,2)
    U_t_scaled(:,i)=U_t(:,i)*sigma_vec_t(i);
end

% compute Ky Fan norm
[V_t,D_t]=eig(X_t'*X_t/n);
[lambda_vec_t,I_t]=sort(diag(D_t),'descend');
kyfan_t=kyfan(lambda_vec_t');

%*************************************
% Applying MCPCA on observed samples

[phi_mat,fun_cell]=MCPCA_sample_mixed_wrapper(X_t,d,q,num_iter,num_init);
X_mc=normalize_matrix(phi_mat);

[U_mc,S_mc,V_mc]=svd(X_mc);
[sigma_vec_mc,I_mc]=sort(diag(S_mc),'descend');
V_mc=V_mc(:,I_mc);
U_mc=U_mc(:,I_mc);
sigma_vec_mc;

% low dimensional embedding using MCPCA
U_mc_scaled=U_mc;
for i=1:size(U_mc,2)
    U_mc_scaled(:,i)=U_mc(:,i)*sigma_vec_mc(i);
end

% compute Ky Fan norm
[V_mc,D_mc]=eig(X_mc'*X_mc/n);
[lambda_vec_mc,I_mc]=sort(diag(D_mc),'descend');
kyfan_mc=kyfan(lambda_vec_mc');

%*************************************
% Comparing quality of low dimensional embeddings

d0=pdist(U00_scaled(:,1:qq));
d_pc=pdist(U_t_scaled(:,1:qq));
d_mc=pdist(U_mc_scaled(:,1:qq));

disp('quality of low dim inference by PCA')
c_pc=corr(d0',d_pc','type','Spearman')

disp('quality of low dim inference by MCPCA')
c_mcpca=corr(d0',d_mc','type','Spearman')


%*************************************
% Comparing KyFan norms

color_mat=[1 0 0
    0 1 0
    0 0 1];

figure
plot(kyfan0,'color',color_mat(1,:))
hold on
plot(kyfan_t,'color',color_mat(2,:))
hold on
plot(kyfan_mc,'color',color_mat(3,:))

xlabel('q')
ylabel('q Kyfan norm')
legend('Latent','PCA','MCPCA')
