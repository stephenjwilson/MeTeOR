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

aucs1=runGSretro_alt('MeTeORgenegene','MeTeORgenegene_2014',AllNetworks, resultdir,bestks{1});
aucs2=runGSretro_alt('MeTeORdiseasegene','MeTeORdiseasegene_2014',AllNetworks, resultdir,bestks{2});
aucs3=runGSretro_alt('MeTeORgenechemical','MeTeORgenechemical_2014',AllNetworks, resultdir,bestks{3});


%% Self-retrospective - MRCOC
if ~exist('MRCOC','var')
    MRCOC = load('MRCOC_2018.mat');
    MRCOC.MRCOC_2015{2} = cellstr(MRCOC.MRCOC_2015{2});
    MRCOC.MRCOC_2018{2} = cellstr(MRCOC.MRCOC_2018{2});
end
aucs4=runGSretro_alt('MRCOC_2018','MRCOC_2015',MRCOC, resultdir,bestks{1});
