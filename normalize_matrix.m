function X=normalize_matrix(X)

% this function makes columsn of X mean zero and norm sqrt(n)

[n,p]=size(X);
%**************
% normalizing the matrix
for i=1:p
    X(:,i)=X(:,i)-mean(X(:,i));
    if norm(X(:,i))~=0
        X(:,i)=X(:,i)./norm(X(:,i))*sqrt(n);
    else
        disp('one of the columns is all zero')
    end
end



