function [ ] = SingleEntityExtraction( meteor,tests,references,k,root,mindata,calc_aucs)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if nargin<6
    mindata=5;
end
if nargin<7
    calc_aucs=1;
end


%% Trim to Gene Only
Names=meteor{1}{2};
Names3=cellfun(@(x) strcmp(x(1:2),'G.'),Names);
Names=Names(Names3);
meteor{1}{2}=Names;
meteor{1}{1}=meteor{1}{1}(Names3,Names3);

%% Map Gene Names
syms={};
inds={};
filename=sprintf('%sgeneMapping.txt',root);
delimiter = '\t';
fileID = fopen(filename,'r');
line_ex = strsplit(fgetl(fileID),delimiter); % header
while ~feof(fileID)
    line_ex = strsplit(fgetl(fileID),delimiter);
    sym=line_ex{3};
    ind=strcat('G.',line_ex{2});
    syms = [syms, sym];
    inds = [inds, ind];
end
fclose(fileID);
[a,b]=ismember(Names,inds);
tmpNames=Names;
tmpNames(a)=syms(b(a));
        
%% Compute Representations
Dinv=sum(meteor{1}{1}).^(-1/2);
norm=diag(Dinv)*meteor{1}{1}*diag(Dinv); %symmetric normalized laplacian L = I - D^(-1/2) * A * D^(-1/2)


%% precompute mappings
gsas = {};
gsa3s = {};
for i=1:length(references)
    GS_=references{i};
    GSnames=GS_{2};
    [a,~]=ismember(GSnames,Names);
    gsas{i}=a;
    [~,a3]=ismember(Names,GSnames(a));
    gsa3s{i}=a3;
end

[ma,~]=ismember(meteor{2}{2},Names);
[~,ma3]=ismember(Names,meteor{2}{2}(ma));

%% Get entity info
aucs = zeros(1,length(Names));
aucsaa = zeros(1,length(Names));
aucscn = zeros(1,length(Names));
aucso = zeros(1,length(Names));
for o=1:length(meteor{2}{2})
    if mod(o,1000)==0
        disp(o);
    end
    entity = meteor{2}{2}{o};
    if ~calc_aucs && ~strcmp(entity,'G.1956')
       continue 
    end
        
    ind=find(strcmp(Names,entity));
    
    SGR=meteor{1}{1}(:,ind);
    if sum(SGR)>2 
        %% Check MsigDB
        GS_=references{2};
        GSnames=GS_{2};
        GS_=GS_{1};
        if sum(strcmp(entity,GSnames))
            ind2=find(strcmp(GSnames,entity));
            GS_=GS_(:,ind2);
            sz=length(Names);
            GS_trim=sparse(1,sz);
            GS_trim(gsa3s{2}>0)=GS_(gsas{2});
            if sum(GS_trim>0)<mindata
                continue
            end
        else
            continue
        end
        
        %% Build Data
        ls={};
        ls{1}=SGR;
        % Add Normalized Vector
        ls{end+1}=norm(:,ind);
        % Add ranks values
        rind=length(ls)+1;
        forwardRank=tiedrank(SGR')';
        meteor_rank = max(forwardRank)-forwardRank+1;
        ls{end+1}= meteor_rank;
        
        % Process MRR
        mrr=MRR(meteor{2}{1},find(strcmp(meteor{2}{2},entity)));
        
        sz=length(Names);
        result_trim=sparse(1,sz);
        
        result_trim(ma3>0)=mrr(ma);
        ls{end+1}=result_trim'; % Trimmed to the right names
        %         ls{end+1}=round(ls{end}); % Rounded
        forwardRank=tiedrank(ls{end}')'; % forwardRank
        ls{end+1}=max(forwardRank)-forwardRank+1; % Rank
        ls{end+1}=ls{rind}-ls{end}; % Relative diff rank
        overall_rank = ls{end};
        
        % Add Test Data
        for o2=1:length(tests)
            SGR=tests{o2}{1}(:,ind);
            if sum(SGR)==0
                continue
            end
            if sum(GS_trim>0)>=mindata && calc_aucs
                tmpSGR = SGR;
                %tmpSGR(ls{1}>0)=0; % zero out MeTeOR
                [~,~,~,AUC] = perfcurve(GS_trim>0,tmpSGR,1);
                if o2==1
                    aucscn(ind)= AUC;
                elseif o2==2
                    aucsaa(ind)= AUC;
                end
            end
            
            ls{end+1}=SGR;
            forwardRank=tiedrank(SGR')'; % forwardRank
            rank_method=max(forwardRank)-forwardRank+1; % Method Rank
            ls{end+1}=meteor_rank-rank_method; % Relative diff rank
            overall_rank = max(overall_rank, meteor_rank-rank_method); % min is max rank
        end
        ls{end+1} = overall_rank;
        %ls{end+1} = 1./overall_rank;
        %ls{end}(ls{1}>0)=0;
        if sum(GS_trim>0)>=mindata && calc_aucs
            tmpSGR = result_trim';
            tmpSGR(ls{1}>0)=0; % zero out MeTeOR
            
            [~,~,~,AUC] = perfcurve(GS_trim>0,tmpSGR,1);
            aucs(ind)= AUC;

            [~,~,~,AUC] = perfcurve(GS_trim>0,1./overall_rank,1);
            aucso(ind)= AUC;
            %[~,b]=ismember(entity,inds);
            %fprintf('id:%s,Sym:%s\n',entity,syms{b});
        end
        
        % Process GSs
        for i=1:length(references)
            GS_=references{i};
            GSnames=GS_{2};
            GS_=GS_{1};
            if sum(strcmp(entity,GSnames))
                ind3=strcmp(GSnames,entity);
                GS_=GS_(:,ind3);
                sz=length(Names);
                GS_trim=sparse(1,sz);
                GS_trim(gsa3s{i}>0)=GS_(gsas{i});
                ls{end+1}=GS_trim'; %#ok<SAGROW>
            else
                ls{end+1}=sparse(length(Names),1);%#ok<SAGROW>
            end
        end

        if strcmp(entity,'G.1956') % only output for EGFR
            %% Write to File
            data=table(tmpNames,ls{:});
            data.Properties.VariableNames = {'Gene','MeTeOR_Weight','MeTeOR_Norm','MeTeOR_Rank','MeTeOR_NMF_Pred','MeTeOR_NMF_Pred_Rank','Diff_Ranking','MeTeOR_CN','MeTeOR_CN_Diff','MeTeOR_AA','MeTeOR_AA_Diff','Max_rank','BIOGRID','MSigDB','STRINGLIT','EVEX'};
            data = sortrows(data,'MeTeOR_Weight','descend');
            writetable(data,sprintf('%svalidation/SingleGene/mrr%s_%s%sv2',root,num2str(k),'MeTeOR',strrep(entity,'.','')),'Delimiter','\t');
        end
    end
end
if calc_aucs
    mask = aucscn>0 | aucsaa>0 | aucs>0 | aucso>0 ;
    data=table(tmpNames(mask),aucs(mask)',aucsaa(mask)',aucscn(mask)',aucso(mask)');
    data.Properties.VariableNames = {'Gene','NMF','AA', 'CN','Combined'};
    data = sortrows(data,'NMF','descend');
    writetable(data,sprintf('%svalidation/SingleGene/%s_%sAUCs',root,num2str(k),'MeTeOR'),'Delimiter','\t');
end
end

