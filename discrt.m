function [y,cen_vec]=discrt(x,B)

% y is a discretized version of x to {1 2 ...B}
% bin sizes are the same

if length(unique(x))>B
    M=max(x);
    m=min(x);
    
    v=linspace(m,M,B+1);
    y=x*0;
    
    cen_vec=zeros(B,1);
    
    for i=1:B
        ind_t=find(x>=v(i) & x<v(i+1));
        y(ind_t)=i;
        cen_vec(i)=(v(i)+v(i+1))/2;
    end
    
    ind_t=find(x==v(B+1));
    y(ind_t)=B;
else
    y=x*0;
    cen_vec=sort(unique(x));
    for i=1:length(cen_vec)
        ind_t=find(x==cen_vec(i));
        y(ind_t)=i;
    end
end


