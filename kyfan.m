function v2=kyfan(v)

n=length(v);
v=sort(v,'descend');
v2(1)=v(1);

for i=2:n
    v2(i)=v2(i-1)+v(i);
end




