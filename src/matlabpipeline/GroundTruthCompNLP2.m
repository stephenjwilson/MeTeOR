function [MCCs,BAUC,BPR,BPR5,BF1,BF15] = GroundTruthCompNLP2( testtitles,comparison,test, AllNetworks,root)
%GROUNDTRUTHCOMP Takes in a list of titles to test against the NLP
%comparison and the MeTeOR test network. 
%   All strings must be titles that are present in All Networks, root is
%   the central location from which results will be stored.
if nargin<4
    load AllNetworks; 
end
if nargin<5
    root='./';
end


if ~exist(root, 'dir')
  mkdir(root);
end

nets=struct();
for i=1:length(testtitles)
    nets.(testtitles{i})=AllNetworks.(testtitles{i});
end


%% Do Comparisons

GStitles = fieldnames(nets);

[titiles1,MCCs1,BAUCs1,BPR1,BPR51,BF11,BF151]=compareGSROC(AllNetworks.(test),GStitles,nets,'MeTeOR',0,root,1); % Normalize  network is a 1 in last option
[titiles2,MCCs2,BAUCs2,BPR2,BPR52,BF12,BF152]=compareGSROC(AllNetworks.(comparison),GStitles,nets,comparison,0,root,0); % Do not normalize references, as they made their normalization or not already


BAUC=zipLists({BAUCs1,BAUCs2});
BPR=zipLists({BPR1,BPR2});
BPR5=zipLists({BPR51,BPR52});
BF1=zipLists({BF11,BF12});
BF15=zipLists({BF151,BF152});
MCCs=zipLists({MCCs1,MCCs2});
Bs={};
for i=1:length(titiles1)
   Bs=[Bs ; titiles1{i}];
   Bs=[Bs ; strcat(titiles2{i},comparison)];
end

title=strcat(root,'Fig2',comparison);
writeMat(strcat(title,'AUCs'),Bs,BAUC);
writeMat(strcat(title,'PRs'),Bs,BPR);
writeMat(strcat(title,'PR5s'),Bs,BPR5);
writeMat(strcat(title,'MCCs'),Bs,MCCs);
end

