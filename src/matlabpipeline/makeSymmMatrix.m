function [ network,e1Names ] = makeSymmMatrix( edgelist,name,root)
%makeSymmMatrix Summary of this function goes here
%   edgelist is a struct that contains at least:
%   edge1,edge2.weight

if nargin<3
    root='./';
end
loc=sprintf('%s%sFAIL.mat',root,name);

try
    load(loc)
catch
    'Failed'
    data=edgelist;
    tmp=data.Entity1;
    data.Entity1=data.Entity2;
    data.Entity2=tmp;
    data=[edgelist;data];
    
    [data] = unique(data,'rows');
    [e1Names,~,e1Indexes] = unique(data.Entity1);
    [e2Names,~,e2Indexes] = unique(data.Entity2);
    
    network=+(sparse(e1Indexes,e2Indexes,data.Confidence,length(e1Names),length(e2Names)));
    n = size(network,1);
    network(1:(n+1):end) = network(1:(n+1):end)/2; % account for the doubling along diagnol
    
    %names=unique([e1names;e2names]);
    
    save(loc,'network','e1Names','-v7.3')
end

end

