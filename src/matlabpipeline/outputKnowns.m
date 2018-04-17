function [] = outputKnowns(AllNetworks,tests,lits,knowns,root)
%UNTITLED Summary of this function goes here
%   tests = a cell array of network names to test
%   lits = a cell array of network names to compare against tests
%   knowns = a cell array of cell arrays with network names of the known
%   information

numKnownTest=zeros(length(tests),1);
numTotalTest=zeros(length(tests),1);
numKnownLit=zeros(length(tests),1);
numTotalLit=zeros(length(tests),1);

for i=1: length(tests)
    [known,total]=getKnown(tests{i},knowns{i},AllNetworks);
    [knownLit,totalLit]=getKnown(lits{i},knowns{i},AllNetworks);
    numKnownTest(i)=known;
    numKnownLit(i)=knownLit;
    numTotalTest(i)=total;
    numTotalLit(i)=totalLit;
end
knownData=table(numKnownTest,numTotalTest,numKnownLit,numTotalLit);
knownData.Properties.VariableNames= {'KnownTest','TotalTest', 'KnownLit', 'TotalLit'};
writetable(knownData,sprintf('%sFig2/KnownData.csv',root))
end

