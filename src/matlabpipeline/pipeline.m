%Pipeline to create data for figures.

%% Check Memory Requirements - Linux Specific. Comment out (but check your memory) if on Windows
% [~,out]=system('vmstat -s -S M | grep "free memory"');
% parts=strsplit(out);
% totalmem=parts{2};
% if str2num(totalmem)<64000
%     error('Total Memory is:%sMB. Recommended Memory is 64GB',totalmem)
% end

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



%% Run Experiments
% Gene Gene Associations
GeneRefs={'MSigDBCurated_top_cp','BIOGRIDCocrystalStructure','BIOGRIDLow','BIOGRIDHigh','BIOGRIDAffinityCaptureWestern','BIOGRIDBiochemicalActivity','BIOGRIDFRET','BIOGRIDReconstitutedComplex'};
GroundTruthCompNLP2(GeneRefs,'EVEX','MeTeORgenegene',AllNetworks,strcat(resultdir,'Fig2/'));
GroundTruthCompNLP2(GeneRefs,'STRING10_textmining','MeTeORgenegene',AllNetworks,strcat(resultdir,'Fig2/'));
% Gene Chemical Associations
ChemRefs=[{'BIOGRIDchem','DGIdb','STITCH5_experimental'} getAllCTD(strcat('../',YamlStruct.general.networkdir))];
GroundTruthCompNLP2(ChemRefs,'STITCH5_textmining','MeTeORgenechemical',AllNetworks,strcat(resultdir,'Fig2/'));
% Gene Disease Associations
DiseaseRefs={'CTD_GD_therapeutic','CTD_GD_marker_mechanism','UNIPROT_fromDisgenet','PSYGENET_fromDisgenet','ORPHANET_fromDisgenet','HPO_fromDisgenet'};
GroundTruthCompNLP2(DiseaseRefs,'BeFree','MeTeORdiseasegene',AllNetworks,strcat(resultdir,'Fig2/'));

%% Get Overlaps
outputKnowns(AllNetworks,{'MeTeORgenegene','MeTeORgenechemical','MeTeORdiseasegene'},{'STRING10_textmining','STITCH5_textmining','BeFree'},{GeneRefs,ChemRefs,DiseaseRefs},resultdir)
getOverlaps(AllNetworks,YamlStruct);
getOverlaps(AllNetworks,YamlStruct,'Lit');

%% Run Sungle Gene

% TP53, KRAS, EGFR, MECP2
Genes={'G.7157','G.3845','G.1956','G.4204'};% Entrez Mapping
GroundTruthComp_SingleGene({'BIOGRIDLow'},'MeTeORgenegene',Genes,AllNetworks,strcat(resultdir,'SingleGene/'));

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


%% Calculate then export EGFR data
% EGFR is 'G.1956' in the mapping
entity='G.1956';
rep=5;
mats={};
Names={};
k=bestks{1};
for i=1:rep
    [ result , Names1]=getNMFresult('MeTeORgenegene',k,AllNetworks,resultdir,3,2,num2str(i));
    mats{i}=result;
    Names{i}=Names1;
end
mrr=MRR(mats,find(strcmp(Names{1},entity)));
references={AllNetworks.BIOGRIDLow,AllNetworks.MSigDBCurated_top_cp,AllNetworks.STRING10_textmining,AllNetworks.EVEX};
Tests={AllNetworks.MeTeORgenegene,{mrr,Names{1}}};
SingleEntityExtraction(Tests,references,entity,k,strcat('../',YamlStruct.general.resultsdir));

exit();
