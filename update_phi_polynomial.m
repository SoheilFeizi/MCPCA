function [phi_mat]=update_phi_polynomial(index,X,d,phi_mat,V)


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

if sum(abs(w))~=0
c_index=polyfit(X(:,index),w,d);
phi_new=polyval(c_index,X(:,index));
else
phi_new=phi_mat(:,index);    
end
%****************************
phi_new=phi_new-mean(phi_new);
if norm(phi_new)~=0
    phi_new=phi_new*sqrt(n)/norm(phi_new);
    phi_mat(:,index)=phi_new;    
else
    % disp('norm phi_new is zero, not updating')
end





