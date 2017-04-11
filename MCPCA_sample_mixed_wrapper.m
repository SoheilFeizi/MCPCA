%function [phi_mat,fun_cell_all]=MCPCA_data_mixed(X_input,d,q,num_iter,num_init)

% INPUTS:
% X_input: input data matrix (n by p), with n samples and p features (discrete or continuous)
% q is the MCPCA rank parameter
% num_iter: number of repeats of coordinate block descent updates
% num_init: number of initializations-1
% d: degrees of freedom for piecewise linear functions for continuous features

% OUTPUTS:
% phi_mat: transformed data matrix (n by p)
% fun_cell: a cell array of size 1 by p, cell i contains tranformed
% alphabet values of vaiable Xi

% Citation:
% Soheil Feizi, David Tse, Maximally Correlated Principal Component Analysis, arXiv:1702.05471

function [phi_mat,fun_cell_all]=MCPCA_data_mixed(X_input,d,q,num_iter,num_init)


[n,p]=size(X_input);
B=d+1;

%****************************
% discretizing data to bins with apx equal number of samples
% use discrt.m to quantize data to bins with equal widths

X_d=X_input*0;
cent_matrix=zeros(B,p);
y_cent_matrix=zeros(B,p);

for i=1:p
    [z1,c1]=discrt_v2(X_input(:,i),B);
    %[z1,c1]=discrt(X_c(:,i),B);
    X_d(:,i)=z1;
    cent_matrix(1:length(c1),i)=c1;
end

[phi_mat,fun_cell_d]=MCPCA_sample_disc_wrapper(X_d,q,num_iter,num_init);

for i=1:p
    map_alphabet=fun_cell_d{i};
    for j=1:B
        ind_t=find(map_alphabet(:,1)==j);
        if length(ind_t)~=0
            y_cent_matrix(j,i)=map_alphabet(ind_t,2);
        else % no mapping value for bin=j
            y_cent_matrix(j,i)=-sqrt(2); % a unique identifier
        end
    end
end

% transform data
phi_mat=X_input*0;
for i=1:p
    ind_t=find(y_cent_matrix(:,i)~=-sqrt(2));
    phi_mat(:,i)=piecewise_lin(X_input(:,i),cent_matrix(ind_t,i),y_cent_matrix(ind_t,i));
end

fun_cell_all=cell(1,p);
for i=1:p
    ind_t=find(y_cent_matrix(:,i)~=-sqrt(2));
    fun_cell_all{1,i}=[cent_matrix(ind_t,i),y_cent_matrix(ind_t,i)];
end


