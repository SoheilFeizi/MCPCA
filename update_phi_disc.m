function [phi_mat,fun_cell]=update_phi_disc_april11(index,X,phi_mat,V,fun_cell)

% updating column index in the phi_mat
% V is a p by q matrix, each column contains an eigenvector

[p,q]=size(V);
n=size(phi_mat,1);

w=zeros(n,1);
for r=1:q
    for i=1:p
        if i~=index
            w=w+V(index,r)*V(i,r)*phi_mat(:,i);
        end
    end
end

%****************************
y_alpha=unique(X(:,index)); % alphabets
M=length(y_alpha); % number of alphabets
w_alpha=zeros(M,1);

for i=1:M
    ind_t=find(X(:,index)==y_alpha(i));
    w_alpha(i,1)=mean(w(ind_t,1));
end

phi_new=zeros(n,1);
for i=1:M
    ind_t=find(X(:,index)==y_alpha(i));
    phi_new(ind_t)=w_alpha(i);
end

%****************************
phi_new=phi_new-mean(phi_new);
if norm(phi_new)~=0
    phi_new=phi_new/norm(phi_new)*sqrt(n);
    phi_mat(:,index)=phi_new;
    fun_cell{index,1}=[y_alpha,w_alpha];
    
else
    % disp('norm phi_new is zero, not updating')
end





