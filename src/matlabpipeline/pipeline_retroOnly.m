%Pipeline to compare naive methods to NMF

%% load yaml library and config
addpath(genpath('yaml'));
yaml_file = '../config.yaml';
YamlStruct = ReadYaml(yaml_file);

%% Set Up Data Folder
 
title=strsplit(YamlStruct.general.title);
resultdir=strcat('../',YamlStruct.general.resultsdir,'validation/');

if ~exist(resultdir, 'dir')
    mkdir(resultdir);
end


%% Load Networks
if exist(strcat('../',YamlStruct.general.networkdir,'AllNetworks.mat'), 'file')==0
    AllNetworks=importNetworks(struct(),YamlStruct);
end

if ~exist('AllNetworks','var')
    load(strcat('../',YamlStruct.general.networkdir,'AllNetworks'));
end


%% CrossValidation - High Memory/Time Requirements!
ks={50,100,200,300,400,500,600};
try
    load(strcat(resultdir,'crossval/bestks.mat'))
catch
    [bestks,bestaucs,MCC1s,MCC10s,MCC50s]=runCrossValidation_symmandnonsymm({'MeTeORgenegene','MeTeORdiseasegene','MeTeORgenechemical'},AllNetworks,strcat(resultdir,'crossval/'),ks);
    save(strcat(resultdir,'crossval/bestks.mat'),'bestks', 'bestaucs','MCC1s','MCC10s','MCC50s');
end

for i=1:3
   r=round(MCC10s(:,i) ,2);
   ind=find(r==max(r),1);
   bestks(i)=ks(ind);
end

%% Self-retrospective - High Memory Requirements if computing NMF predictions from scratch

[sp1,BAUCs1,BPR1,BPR51,BF11,BF151]=runGSretro('MeTeORgenegene','MeTeORgenegene_2014','MeTeORgenegene','MeTeORgenegene_2014',AllNetworks, resultdir,bestks{1});
[sp2,BAUCs2,BPR2,BPR52,BF12,BF152]=runGSretro('MeTeORdiseasegene','MeTeORdiseasegene_2014','MeTeORdiseasegene','MeTeORdiseasegene_2014',AllNetworks, resultdir,bestks{2});
[sp3,BAUCs3,BPR3,BPR53,BF13,BF153]=runGSretro('MeTeORgenechemical','MeTeORgenechemical_2014','MeTeORgenechemical','MeTeORgenechemical_2014',AllNetworks, resultdir,bestks{3});


aps={sp1{:},sp2{:},sp3{:}};
BAUC=[BAUCs1,BAUCs2,BAUCs3];
BPR=[BPR1,BPR2,BPR3];
BPR5=[BPR51,BPR52,BPR53];
BF1=[BF11,BF12,BF13];
BF15=[BF151,BF152,BF153];

title='selfretro';
boxroot=strcat(resultdir,'/',title,'/BoxPlot/');
if ~exist(boxroot, 'dir')
    mkdir(boxroot);
end

plotBoxPlot(BAUC,aps,strcat(title,'AUC'),boxroot);
plotBoxPlot(BPR,aps,strcat(title,'PR1'),boxroot);
plotBoxPlot(BPR5,aps,strcat(title,'PR5'),boxroot);
plotBoxPlot(BF1,aps,strcat(title,'F11'),boxroot);
plotBoxPlot(BF15,aps,strcat(title,'F15'),boxroot);

fl=strcat(resultdir,'/',title,'/');
writeMatnoMod(strcat(fl,'AUCs'),aps,BAUC);
writeMatnoMod(strcat(fl,'PRs'),aps,BPR);
writeMatnoMod(strcat(fl,'PR5s'),aps,BPR5);
