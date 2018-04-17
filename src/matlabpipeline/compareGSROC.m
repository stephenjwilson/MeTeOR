function [ titles,MCCs,BAUCs,BPR,BPR5,BF1,BF15] = compareGSROC(test,GStitles,GSnet,title,initial,root, norm)
%compareGSROC Does a comparative analysis of multiple networks to Gold
% Standards.
%   Test is cell array with the network and labels
%   GStitles is a cell array of strings that match to networks in GSnet
%   Title is the name of the experiment
%   initial is a binary focusing 
if nargin<7
    norm=1;
end
if nargin<6
    root='./';
end
if nargin<5
    initial=0;
else
    if initial
        title=sprintf('%s_initial',title)
    end
end

try
    title
    load(sprintf('%s.mat_FAIL',title))
    disp('Loaded From File')
catch
    %Initialize vectors

    MCCs=zeros(1,length(GStitles));
    BAUCs=zeros(1000,length(GStitles));
    BPR=zeros(1000,length(GStitles));
    BPR5=zeros(1000,length(GStitles));
    BF1=zeros(1000,length(GStitles));
    BF15=zeros(1000,length(GStitles));
    
    Testnetwork=test{1};
    if norm
        Dinv=sum(Testnetwork).^(-1/2);
        Testnetwork=diag(Dinv)*Testnetwork*diag(Dinv); %symmetric normalized laplacian L = I - D^(-1/2) * A * D^(-1/2)
    end
    Testnames=test{2};
    [~,h]=size(Testnames);

    if h==2
        Testnames1=Testnames{1};
        Testnames2=Testnames{2};
    else
        Testnames1=Testnames;
        Testnames2=Testnames;
    end    

    if min(min(Testnetwork))~=0
       Testnetwork(Testnetwork~=0)=Testnetwork(Testnetwork~=0)+abs(min(min(Testnetwork)));
    end    
    
    for i=1:length(GStitles)
        netname=GStitles{i};
        GS=GSnet.(netname){1};
        GSnames=GSnet.(netname){2};
        [~,h]=size(GSnames);
        
        if h==2
            GSnames1=GSnames{1};
            GSnames2=GSnames{2};
        else
            GSnames1=GSnames;
            GSnames2=GSnames;
        end

        %% Get unique genes, elminiate empty strings
        [Testnetwork,Testnames1,Testnames2]=eliminateDuplicateNames(Testnetwork,Testnames1,Testnames2,1);
        [GS,GSnames1,GSnames2]=eliminateDuplicateNames(GS,GSnames1,GSnames2,1);
        
        %% Figure out mapping to a common set of proteins
        [Test_to_common_mapping1]=ismember(Testnames1,GSnames1);
        
        commonNames1=Testnames1(Test_to_common_mapping1);
        [GS_to_common_mapping1]=ismember(GSnames1,commonNames1);
        
        [Test_to_common_mapping2]=ismember(Testnames2,GSnames2);
        commonNames2=Testnames2(Test_to_common_mapping2);
        [GS_to_common_mapping2]=ismember(GSnames2,commonNames2);
        
        disp(sprintf('Found Mappings: ax1 %s/%s,ax2 %s/%s',num2str(length(commonNames1)),num2str(length(Testnames1)),num2str(length(commonNames2)),num2str(length(Testnames2))));
        
        
        %% MCCs - NOT WORKING
        %pred=Testnetwork(Test_to_common_mapping1,Test_to_common_mapping2);
        %labels=GS(GS_to_common_mapping1,GS_to_common_mapping2);
        %MCC=calcMCC(labels,pred);
        MCCs(i)=0;
        %% Bootstrap
        [aucs,PRs,PR5s,F1s,F15s ]=bootstrapComparison(Testnetwork(Test_to_common_mapping1,Test_to_common_mapping2),GS(GS_to_common_mapping1,GS_to_common_mapping2),strcat(root,'/BoxPlot/'),strcat(title,netname));
        BAUCs(1:nnz(aucs),i)=aucs;
        BPR(1:nnz(aucs),i)=PRs;
        BPR5(1:nnz(aucs),i)=PR5s;
        BF1(1:nnz(aucs),i)=F1s;
        BF15(1:nnz(aucs),i)=F15s;


    end    
    %stats
    %save(sprintf('%s.mat',title),'-v7.3')
end

%Determine Graph Labels

f2=@(name) sprintf('%s%s', title,name);
titles=cellfun(f2, GStitles', 'UniformOutput',false);

% %% MetaBoxplot
% boxroot=strcat(root,'/BoxPlot/');
% plotBoxPlot(BAUCs,s2,strcat(title,'AUC'),boxroot,1:length(s2));
% plotBoxPlot(BPR,s2,strcat(title,'PR1'),boxroot,1:length(s2));
% plotBoxPlot(BPR5,s2,strcat(title,'PR5'),boxroot,1:length(s2));
% plotBoxPlot(BF1,s2,strcat(title,'F11'),boxroot,1:length(s2));
% plotBoxPlot(BF15,s2,strcat(title,'F15'),boxroot,1:length(s2));
% 
% disp('Plotted')

end

