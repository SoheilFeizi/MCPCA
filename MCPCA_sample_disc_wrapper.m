%function [phi_mat,fun_cell]=MCPCA_sample_disc_wrapper(X_input,q,num_iter,num_init)

% INPUTS:
% X_input: input data matrix (n by p), with n samples and p features
% q is the MCPCA rank parameter
% num_iter: number of repeats of coordinate block descent updates
% num_init: number of initializations-1

% OUTPUTS:
% phi_mat: transformed data matrix (n by p)
% fun_cell: a cell array of size 1 by p, cell i contains tranformed
% alphabet values of vaiable Xi

% Citation:
% Soheil Feizi, David Tse, Maximally Correlated Principal Component Analysis, arXiv:1702.05471


function [phi_mat,fun_cell]=MCPCA_sample_disc_wrapper(X_input,q,num_iter,num_init)


%*********************************
% Labeling alphabets with {1,2,...,B_i}

X0=X_input*0;
[n,p]=size(X0);
map_alphabet0=cell(1,p);
for i=1:p
    alpha0=unique(X_input(:,i));
    alpha1=(1:length(alpha0))';
    map_alphabet0{1,i}=[alpha0,alpha1];
    for j=1:length(alpha0)
        ind_t=find(X_input(:,i)==alpha0(j));
        X0(ind_t,i)=j;
    end
end

%*********************************
% Applying MCPCA with identity initialization

disp('MCPCA: iteration 0')
[obj_vec_test,phi_mat_test,~]=MCPCA_sample_disc(X0,q,num_iter,0);

%*********************************
% extracting functions

fun_cell=cell(p,1);
for ii=1:p
    alpha0=unique(X_input(:,ii));
    alpha1=alpha0*0;
    for jj=1:length(alpha0)
        ind_t=find(X_input(:,ii)==alpha0(jj));
        if length(unique(phi_mat_test(ind_t,ii)))==1
            alpha1(jj)=unique(phi_mat_test(ind_t,ii));
        end
    end
    fun_cell{ii,1}=[alpha0,alpha1];
end
%*********************************
% saving results for init 0

total_res=cell(3,num_init+1);
kyfan_vec=zeros(1,num_init+1);

kyfan_vec(1)=obj_vec_test(end);
total_res{1,1}=obj_vec_test;
total_res{2,1}=phi_mat_test;
total_res{3,1}=fun_cell;


%*********************************
% random initialization of MCPCA

for i=1:num_init
    disp(['MCPCA: iteration ',num2str(i)])
    [obj_vec_test,phi_mat_test,~]=MCPCA_sample_disc(X0,q,num_iter,1);
    
    % extracting functions
    fun_cell=cell(p,1);
    for ii=1:p
        alpha0=unique(X_input(:,ii));
        alpha1=alpha0*0;
        for jj=1:length(alpha0)
            ind_t=find(X_input(:,ii)==alpha0(jj));
            if length(unique(phi_mat_test(ind_t,ii)))==1
                alpha1(jj)=unique(phi_mat_test(ind_t,ii));
            end
        end
        fun_cell{ii,1}=[alpha0,alpha1];
    end
    
    % saving results for init i
    kyfan_vec(i+1)=obj_vec_test(end);
    total_res{1,i+1}=obj_vec_test;
    total_res{2,i+1}=phi_mat_test;
    total_res{3,i+1}=fun_cell;
end

%*********************************
% selecting best initilization
[~,ind_m]=max(kyfan_vec);
phi_mat=total_res{2,ind_m};
fun_cell=total_res{3,ind_m};

%*********************************
% making positive orientation
[n,p]=size(X_input);
for i=1:p
    temp=corrcoef(X_input(:,i),phi_mat(:,i));
    if temp(1,2)<0
        phi_mat(:,i)=-phi_mat(:,i);
        fun_temp=fun_cell{i,1};
        fun_temp(:,2)=-fun_temp(:,2);
        fun_cell{i,1}=fun_temp;
    end
end