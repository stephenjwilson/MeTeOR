function [ mrr ] = MRR( mats,index, axis)
%MRR Summary of this function goes here
%   Detailed explanation
if nargin<3
    axis=0;
end
m=length(mats);
[h,w]=size(mats{1});

if axis
    data=vertcat(mats{:})';
    n=h;
else
    data=horzcat(mats{:});
    n=w;
end

vec=data(index,:);
mat=reshape(vec,n,m)';

forwardRank=tiedrank(mat')';
rank=repmat(max(forwardRank,[],2),1,w)-forwardRank+1;

mrr=sum(1./rank)/w;
mrr(mrr==min(mrr))=0;
end

