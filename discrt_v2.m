function [y,cen_vec]=discrt_v2(x,B)

% y is a discretized version of x to {1 2 ...B}
% each bin has apx similar number of samles

x_m=mode(x);
n=length(x);

if length(unique(x))>B && sum(x==x_m)<floor(n/B)
    %continuous case
    [x_s,I]=sort(x);
    n=length(x_s);
    M=floor(n/B); % number of points in each bin
    y_s=x_s*0;
    cen_vec=zeros(B,1);
    y_s(1:M)=1;
    cen_vec(1)=mean(x_s(1:M));
    
    for i=2:B-1
        if mean(x_s((i-1)*M+1:i*M))~=cen_vec(i-1)
            cen_vec(i)=mean(x_s((i-1)*M+1:i*M));
            y_s((i-1)*M+1:i*M)=i;
        else
            disp('use discrt.m instead of discrete_v2')
            return;
        end
    end
    y_s((B-1)*M+1:end)=B;
    cen_vec(B)=mean(x_s((B-1)*M+1:end));
    
    y=y_s*0;
    y(I)=y_s;
    
elseif length(unique(x))>B && sum(x==x_m)>=floor(n/B)
    % mixed disc/con
    
    x_uni=sort(unique(x));
    n=length(x);
    M=floor(n/B);
    
    f_x_uni=x_uni*0;
    for i=1:length(x_uni)
        f_x_uni(i)=sum(x==x_uni(i));
    end
    
    
    f_ind=f_x_uni*0;
    f_ind(1)=1;
    for i=2:length(f_ind)
        
        ind_t=find(f_ind==f_ind(i-1));
        %sum(f_x_uni(ind_t));% all elements with previous index
        if sum(f_x_uni(ind_t))>=M
            f_ind(i)=f_ind(i-1)+1;% new index
        else
            f_ind(i)=f_ind(i-1);
        end
    end
    
    y=x*0;
    cen_vec=zeros(1,B);
    
    for i=1:B
        ind_t=find(f_ind==i);
        if length(ind_t)~=0
            cen_vec(i)=mean(x_uni(ind_t).*f_x_uni(ind_t));
            
            for j=1:length(ind_t)
                ind_tt=find(x==x_uni(ind_t(j)));
                y(ind_tt)=i;
            end
        end
    end
       
elseif length(unique(x))<=B
    y=x*0;
    cen_vec=sort(unique(x));
    for i=1:length(cen_vec)
        ind_t=find(x==cen_vec(i));
        y(ind_t)=i;
    end
end



