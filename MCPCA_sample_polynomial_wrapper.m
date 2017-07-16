%function [phi_mat,fun_cell_all]=MCPCA_sample_polynomial_wrapper(X_input,d,q,num_iter,num_init)

% INPUTS:
% X_input: input data matrix (n by p), with n samples and p features (discrete or continuous)
% q is the MCPCA rank parameter
% num_iter: number of repeats of coordinate block descent updates
% num_init: number of initializations-1
% d: maximum degree of each polynomial for continuous features


% OUTPUTS:
% phi_mat: transformed data matrix (n by p)
% fun_cell: a cell array of size 1 by p, cell i contains polynomial coefs of vaiable Xi

% Citation:
% Soheil Feizi, David Tse, Maximally Correlated Principal Component Analysis, arXiv:1702.05471

function [phi_mat,fun_cell]=MCPCA_sample_polynomial_wrapper(X_input,d,q,num_iter,num_init)

%*********************************
% Applying MCPCA with identity initialization

X0=normalize_matrix(X_input);

disp('MCPCA: initialization 0')
[obj_vec_test,phi_mat_test]=MCPCA_sample_polynomial(X0,d,q,num_iter,0);

%*********************************
% saving results for init 0

total_res=cell(2,num_init+1);
kyfan_vec=zeros(1,num_init+1);

kyfan_vec(1)=obj_vec_test(end);
total_res{1,1}=obj_vec_test;
total_res{2,1}=phi_mat_test;


%*********************************
% random initialization of MCPCA

for i=1:num_init
    disp(['MCPCA: initialization ',num2str(i)])
    [obj_vec_test,phi_mat_test]=MCPCA_sample_polynomial(X0,d,q,num_iter,1);
    
    % saving results for init i
    kyfan_vec(i+1)=obj_vec_test(end);
    total_res{1,i+1}=obj_vec_test;
    total_res{2,i+1}=phi_mat_test;
end

%*********************************
% selecting best initilization
[~,ind_m]=max(kyfan_vec);
obj_vec=total_res{1,ind_m};
phi_mat=total_res{2,ind_m};

%*********************************
[n,p]=size(X_input);
for i=1:p
    temp=corrcoef(X_input(:,i),phi_mat(:,i));
    if temp(1,2)<0
        phi_mat(:,i)=-phi_mat(:,i);
    end
end

%****************************
phi_mat=normalize_matrix(phi_mat);
fun_cell=cell(p,1);

for index=1:p
    cc=polyfit(X_input(:,index),phi_mat(:,index),d);
    fun_cell{index,1}=cc;
end





