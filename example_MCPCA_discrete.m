% example script
% MCPCA on discrete data

% Citation:
% Soheil Feizi, David Tse, Maximally Correlated Principal Component Analysis, arXiv:1702.05471


clc
close all
clear all

%*************************************
% Parameters of generaing data

n=1000;% number of samples
p=50; %number of variables
B=10; %size of variable alphabets

K=4; % number of clusters in feature covariance matrix
p_in=0.75; % corvariance matrix parameters
p_across=0; % corvariance matrix parameters

%*************************************
%MCPCA parameters

q=4;
num_init=10; % number of initializations
num_iter=10; % number of repeats of coordinate block descent updates

%*************************************
% Generate latent data samples

ind_cl=ceil(linspace(1,p,K+1));
ind_cl_true=zeros(p,1);
clust=cell(1,K);
for i=1:K
    if i~=K
        clust{1,i}=ind_cl(i):(ind_cl(i+1)-1);
        ind_cl_true(ind_cl(i):(ind_cl(i+1)-1))=i;
    else
        clust{1,i}=ind_cl(i):(ind_cl(i+1));
        ind_cl_true(ind_cl(i):(ind_cl(i+1)))=i;
    end
end

cluster_matrix_true=zeros(p);
for i=1:K
    ind_t=find(ind_cl_true==i);
    cluster_matrix_true(ind_t,ind_t)=1;
end

cov_mat=make_sbm(p,clust,p_in,p_across);
[V,D]=eig(cov_mat);
[lambda_vec,I]=sort(diag(D),'descend');
V=V(:,I);

cov_mat_sdp=cov_mat*0;
for i=1:p
    if lambda_vec(i)>0
        cov_mat_sdp=cov_mat_sdp+lambda_vec(i)*V(:,i)*V(:,i)';
    end
end

X0=mvnrnd(zeros(1,p),cov_mat_sdp,n);
X0=normalize_matrix(X0);

% discretizing data
X_d=X0*0;
for i=1:p
    [z,~]=discrt(X0(:,i),B);
    X_d(:,i)=z;
end

X_d=normalize_matrix(X_d);
[V_d,D_d]=eig(X_d'*X_d/n);

[lambda_vec_d,I_d]=sort(diag(D_d),'descend');
V_d=V_d(:,I_d);
lambda_mat_all=lambda_vec_d';

figure
imagesc(abs(X_d'*X_d/n))
title('latent low rank abs covariance')
colorbar

%*************************************
% Generate observed data samples
% randomly transforming data

X=X_d*0;
for i=1:p
    alphabet=unique(X_d(:,i));
    alphabet2=randperm(length(alphabet));
    for j=1:length(alphabet)
        ind_t=find(X_d(:,i)==alphabet(j));
        X(ind_t,i)=alphabet2(j);
    end
end

% making variables mean zero unit variance
X_t=normalize_matrix(X);

[V_t,D_t]=eig(X_t'*X_t/n);
[lambda_vec_t,I_t]=sort(diag(D_t),'descend');
V_t=V_t(:,I_t);
lambda_mat_all=[lambda_mat_all;lambda_vec_t'];

figure
imagesc(abs(X_t'*X_t/n))
title('observed abs covariance')
colorbar

%*************************************
% Applying MCPCA on observed samples

[phi_mat,fun_cell]=MCPCA_sample_disc_wrapper(X,q,num_iter,num_init);

phi_mat=normalize_matrix(phi_mat);
[V_MCPCA,D_MCPCA]=eig(phi_mat'*phi_mat/n);

[lambda_vec_MCPCA,I_MCPCA]=sort(diag(D_MCPCA),'descend');
V_MCPCA=V_MCPCA(:,I_MCPCA);
lambda_mat_all=[lambda_mat_all;lambda_vec_MCPCA'];

figure
imagesc(abs(phi_mat'*phi_mat/n))
title(['MCPCA abs covariance-',num2str(q)])
colorbar

%*************************************
% Comparing KyFan norms

color_mat=[1 0 0
    0 1 0
    0 0 1];

figure
for i=1:3
    v_t=lambda_mat_all(i,:);
    v_t2=kyfan(v_t);
    plot(v_t2,'color',color_mat(i,:))
    hold on
end
xlabel('q')
ylabel('q Kyfan norm')
legend('Latent','PCA','MCPCA')

