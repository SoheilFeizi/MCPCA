function y=piecewise_lin(x,x_c,y_c)

% x_c is the vector of length B
% y_c is a vector of length B, corersponding to values of x_c
% x is a vector of length n.
% y: values at x using a piecewise linear fit of (x_c,y_c)

y=x*0;
B=length(x_c);
if B<2
    return;
    disp('min alhabet size is 2')
end

% sort x_c
[x_c,I]=sort(x_c);
y_c=y_c(I);

% the lower tail
ind0=find(x<=x_c(2));
m0=(y_c(2)-y_c(1))/(x_c(2)-x_c(1));
y0=y_c(1)-m0*x_c(1);
y(ind0)=m0*x(ind0)+y0;

% the upper tail
ind0=find(x>=x_c(B-1));
m0=(y_c(B)-y_c(B-1))/(x_c(B)-x_c(B-1));
y0=y_c(B)-m0*x_c(B);
y(ind0)=m0*x(ind0)+y0;

if B>=4
    for i=2:B-2
        ind0=find(x<=x_c(i+1) & x>=x_c(i));
        m0=(y_c(i+1)-y_c(i))/(x_c(i+1)-x_c(i));
        y0=y_c(i)-m0*x_c(i);
        y(ind0)=m0*x(ind0)+y0;
        
    end
end


