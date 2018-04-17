function [numKnown,totalNum] = getKnown(test,known,AllNetworks)
%UNTITLED Summary of this function goes here
%   test = a title for the test network
%   known = a cell array of titles for the knowns to compare against
%   AllNetworks = a struct that contains all network information

net=AllNetworks.(test){1};
names=AllNetworks.(test){2};

totalNum=nnz(net);
[w,h]=size(net);
knownMat=sparse(w,h);

for i=1:length(known)
    knownTitle=known{i}; 
    netKnown=AllNetworks.(knownTitle){1};
    namesKnown=AllNetworks.(knownTitle){2};
    %% Figure out mapping to a common set of proteins
    [a,b]=ismember(namesKnown,names);
    %% Get data from known into test-like matrix
    tmp=sparse(w,h);
    tmp(b(b>0),b(b>0))=netKnown(a,a);
    %% Intersection of known with test and then Union with previous known data 
    knownMat=or(knownMat,and(tmp,net));
end

numKnown=nnz(knownMat);
end

