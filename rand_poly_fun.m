function [y,c_f]=rand_poly_fun(x,d)

% x is a vector of length n
% we apply a random polynomial with degree d to x
% c_f is the coef vector of length d+1

n=length(x);
c=randn(d+1,1);
c=c/norm(c);

y=x*0;
for i=0:d
    y=y+c(i+1)*hermite(i,x);
end

y=y-mean(y);
y=y/norm(y)*sqrt(n);


c_f=polyfit(x,y,d);



