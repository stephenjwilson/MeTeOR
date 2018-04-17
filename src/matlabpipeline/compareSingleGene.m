function [ X,Y,s1,XP,YP,s2] = compareSingleGene(test,GStitles,GSnet,focus, title,initial,root, norm)
%compareGSROC Does a comparative analysis of multiple networks to a Gold
% Standard. Title is the name of the GS, testtitles is a cell array of the
% names of the methods being tested

if nargin<8
    norm=1;
end
if nargin<7
    root='./';
end
if nargin<6
    initial=0;
else
    if initial
        title=sprintf('%s_initial',title)
    end
end
title=sprintf('%s_%s',title,focus)
try
    title
    load(sprintf('%s.mat_FAIL',title))
    disp('Loaded From File')
catch
    %Initialize vectors
    X=zeros(1000,length(GStitles));
    Y=zeros(1000,length(GStitles));
    XP=zeros(1000,length(GStitles));
    YP=zeros(1000,length(GStitles));
    AUC=zeros(1,length(GStitles));
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
%     if ~isempty(strfind(title,'EVEX'))
% %         Testnetwork=Testnetwork*-1;
% %         m=min(min(Testnetwork));
% %         ind=find(Testnetwork);
% %         Testnetwork(ind)=Testnetwork(ind)-m;
%         Testnetwork=abs(Testnetwork);
%     end
    if min(min(Testnetwork))~=0
       Testnetwork(Testnetwork~=0)=Testnetwork(Testnetwork~=0)+abs(min(min(Testnetwork)));
    end    
    
    for i=1:length(GStitles)
        netname=GStitles{i};
        GS=GSnet.(netname){1};
        GSnames=GSnet.(netname){2};
        [w,h]=size(GSnames);
        if min(min(GS))~=0 % This thresholds Reactome 
           GS(GS<2)=0;
        end   
        
        if h==2
            GSnames1=GSnames{1};
            GSnames2=GSnames{2};
        else
            GSnames1=GSnames;
            GSnames2=GSnames;
        end
        
%         if strfind(netname,'STRING')
%             GS=GS>900;
%         end
        
        %% Get unique genes, elminiate genes that don't have HGNC ids (empty strings)
%         Testnames1=unique(Testnames1);
%         Testnames1(find(strcmp('',Testnames1)))=[]; %#ok<FNDSB>
%         Testnames2=unique(Testnames2);
%         Testnames2(find(strcmp('',Testnames2)))=[]; %#ok<FNDSB>
%         GSnames1=unique(GSnames1);
%         GSnames1(find(strcmp('',GSnames1)))=[]; %#ok<FNDSB>
%         GSnames2=unique(GSnames2);
%         GSnames2(find(strcmp('',GSnames2)))=[]; %#ok<FNDSB>

        %[Testnetwork,Testnames1,Testnames2]=eliminateDuplicateNames(Testnetwork,{focus},Testnames2,1);
        %[GS,GSnames1,GSnames2]=eliminateDuplicateNames(GS,{focus},GSnames2,1);
        
        %% Figure out mapping to a common set of proteins
        [Test_to_common_mapping1]=ismember(Testnames1,GSnames1);
        
        commonNames1=Testnames1(Test_to_common_mapping1);
        [GS_to_common_mapping1]=ismember(GSnames1,commonNames1);
        
        [Test_to_common_mapping2]=ismember(Testnames2,GSnames2);
        commonNames2=Testnames2(Test_to_common_mapping2);
        [GS_to_common_mapping2]=ismember(GSnames2,commonNames2);
        
        Testnetwork_trim=Testnetwork(strcmp(Testnames1,focus),Test_to_common_mapping1);
        GS_trim=GS(strcmp(GSnames1,focus),GS_to_common_mapping1);
        nnz(and(GS_trim>0,Testnetwork_trim>0))
        %disp(sprintf('Found Mappings: ax1 %s/%s,ax2 %s/%s',num2str(length(commonNames1)),num2str(length(Testnames1)),num2str(length(commonNames2)),num2str(length(Testnames2))));
        
        %% Find positive common edges from test network, flatten matrix
        
        positives=and(Testnetwork_trim>0,GS_trim>0); % If the trimmed GS is nonzero and trimmed test is nonzero, it is a positive
        positives=positives(Testnetwork_trim>0);
        Testnetwork_trim=Testnetwork_trim(Testnetwork_trim>0);
        %% Get ROC curve
        [x_ROC,y_ROC,~,auc_ROC]=perfcurve(positives,Testnetwork_trim,true);
        [TP,FP,~,~]=perfcurve(positives,Testnetwork_trim, true, 'XCrit','TP','YCrit','FP');
        fprintf('TP:%d,FP:%d,sum:%d,nnz:%d\n',TP(end),FP(end), sum([TP(end),FP(end)]),nnz(Testnetwork_trim) );
        
%         [TN,FN,~,~]=perfcurve(positives_sorted_trim,Testnetwork_sorted_trim, true, 'XCrit','TN','YCrit','FN');
%         fprintf('TN:%d,FN:%d,sum:%d,nnz:%d\n',TN(1),FN(1), sum([TN(1),FN(1)]),nnz(Testnetwork2) );
        
        disp(auc_ROC)

        
        %% Reduce dimensionality if needed
        if length(x_ROC)>1000
            x_ROC=reduceVector(x_ROC,1000);
            y_ROC=reduceVector(y_ROC,1000);
        else
            x_ROC(end:1000)=1;
            y_ROC(end:1000)=1;
        end

        
        %% Precision recall curves
        [x_precRecall,y_precRecall,~,~]=perfcurve(positives,Testnetwork_trim, true, 'XCrit','TPR','YCrit','PPV');

        if length(x_precRecall)>1000
            x_precRecall=reduceVector(x_precRecall,1000);
            y_precRecall=reduceVector(y_precRecall,1000);
        else
            x_precRecall(end:1000)=x_precRecall(end);
            y_precRecall(end:1000)=y_precRecall(end);
        end

        
        r=0.1;
        ind=find(x_precRecall>r,1);
        p=y_precRecall(ind);
        r=x_precRecall(ind);
        f1=2*(p*r)/(p+r);
        r=0.05;
        ind=find(x_precRecall>r,1);
        p=y_precRecall(ind);
        f105=2*(p*r)/(p+r);
        disp(sprintf('F1 at 0.1:%f,F1 at 0.5:%f',f1,f105));
        %disp('Got PR')
        
        %% Store all information
        X(:,i)=x_ROC;
        Y(:,i)=y_ROC;
        XP(:,i)=x_precRecall;
        YP(:,i)=y_precRecall;
        AUC(i)=auc_ROC;
        

    end    
    %stats
    %save(sprintf('%s.mat',title),'-v7.3')
end

%Determine Graph Labels


f1=@(name,a) sprintf('%s %s AUC=%0.2f',name, focus,a);
f2=@(name) sprintf('%s %s', name, focus);

s1=cellfun(f1, GStitles', num2cell(AUC),'UniformOutput',false);
s2=cellfun(f2, GStitles','UniformOutput',false);

%plot curves
plot_perfcurve(X,Y,s1,0,title,root);
plot_perfcurve(XP,YP,s2,1,title,root);


disp('Plotted')

end

