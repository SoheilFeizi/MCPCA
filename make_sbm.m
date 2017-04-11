function [ adj ] = make_sbm( N,clust,p,q )

% make_sbm draws a graph from the stochastic block model with given
% parameters
% N: graph size;
% clust: a cell array specifying K communities;
% p: within-cluster density,
% q: scalar outside-of-cluster density.

K=numel(clust);

% error handling
for k=1:K
    if size(clust{k},1)~=1
        error('specify clusters as row vectors in a cell array')
    end
end

v=sort([clust{:}]);
if numel(v)==N
    if any(v~=1:N)
        error('Each node i=1,...,N should be in exactly one community')
    end
else
    error('Each node i=1,...,N should be in exactly one community')
end


% vector of within-cluster densities
if numel(p)==1
    p=ones(numel(clust),1)*p;
end

% generate adjacency matrix
adj=rand(N)<q;
for k=1:numel(clust)
    adj(clust{k},clust{k})=rand(numel(clust{k}))<p(k);
end
adj=triu(adj);
adj=adj+adj';
adj=adj-diag(diag(adj));

end

