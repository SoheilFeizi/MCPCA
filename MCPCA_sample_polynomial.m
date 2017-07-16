function [obj_vec,phi_mat]=MCPCA_sample_polynomial(X,d,q,num_iter,control_init)

% INPUTS:
% X_input: input data matrix (n by p), with n samples and p features
% q is the MCPCA rank parameter
% num_iter: number of repeats of coordinate block descent updates
% control_init=0, it starts from phi(X)=X
% control_init=1 it starts from random phi(X)


% OUTPUTS:
% phi_mat: transformed data matrix (n by p)
% obj_vec: q Ky Fan norm

% Citation:
% Soheil Feizi, David Tse, Maximally Correlated Principal Component Analysis, arXiv:1702.05471

%*********************************
% initialization
[n,p]=size(X);

if control_init==0
    phi_mat=X;
else
    phi_mat=X*0;
    for i=1:p
        phi_mat(:,i)=rand_poly_fun(X(:,i),d);
    end
end

phi_mat=normalize_matrix(phi_mat);

Z=phi_mat'*phi_mat/n;
[V,D]=eigs(Z,q);
lam_vec=diag(D);

obj_vec=zeros(1,num_iter+1);
obj_vec(1,1)=sum(lam_vec);


for iter=1:num_iter
    %updating phi_i
    for i=1:p
        [phi_mat]=update_phi_polynomial(i,X,d,phi_mat,V);
    end
    
    % updating V variables
    Z=phi_mat'*phi_mat/n;
    [V,D]=eigs(Z,q);
    lam_vec=diag(D);
    
    % objective value
    obj_vec(1,iter+1)=sum(lam_vec);
end




