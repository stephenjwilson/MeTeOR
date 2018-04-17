function [ x,y,s,xp,yp,s2] = GroundTruthComp_SingleGene( testtitles,test,Entities,AllNetworks,root)
%GROUNDTRUTHCOMP Creates an ROC curve based on a single-gene's rankings
%   Detailed explanation goes here
if ~exist(root, 'dir')
  mkdir(root);
end
x=[];
y=[];
s={};
xp=[];
yp=[];
s2={};

nets=struct();
for i=1:length(testtitles)
    nets.(testtitles{i})=AllNetworks.(testtitles{i});
end

GStitles = fieldnames(nets);
for i=1:length(Entities)
    [X,Y,S,XP,YP,S2]=compareSingleGene(AllNetworks.(test),GStitles,nets,Entities{i},test,0,root,0);
    x=[x X];
    y=[y Y];
    s=[s S];
    xp=[xp XP];
    yp=[yp YP];
    s2=[s2 S2];
end

title='SingleGene';
%plot curves
plot_perfcurve(x,y,s,0,title,root);
plot_perfcurve(xp,yp,s2,1,title,root);

end

