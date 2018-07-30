function [ ] = SingleEntityExtractionNaive(Tests,References,entity,root)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Trim to Gene Only
Names=Tests{1}{2};
Names3=cellfun(@(x) strcmp(x(1:2),'G.'),Names);
Names=Names(Names3);
Tests{1}{2}=Names;
Tests{1}{1}=Tests{1}{1}(Names3,Names3);

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

%% Get entity info

ind=find(strcmp(Names,entity));

SGR=Tests{1}{1}(:,ind);
if sum(SGR)>0
    ls={};
    ls{1}=SGR;
    
    % Process GSs
    for i=1:length(References)
        GS_=References{i};
        GSnames=GS_{2};
        GS_=GS_{1};
        if sum(strcmp(entity,GSnames))
            ind=find(strcmp(GSnames,entity));
            GS_=GS_(:,ind);
            [a,~]=ismember(GSnames,Names);
            sz=length(Names);
            GS_trim=sparse(1,sz);
            [~,a3]=ismember(Names,GSnames(a));
            GS_trim(a3>0)=GS_(a);
            ls{end+1}=GS_trim'; %#ok<SAGROW>
        else
            ls{end+1}=sparse(length(Names),1);%#ok<SAGROW>
        end
    end
    
    %data=table(Names,SGMeTeOR,SGMeTeOR_r,SGMeTeOR_fact,SGMeTeOR_fact_r,ls{1},ls{2},ls{3});
    %data.Properties.VariableNames = {'Gene','Weight','WeightRank','WeightFact','WeightFactRank','STRING','EVEX','Reactome'};
    
    %% Map to symbols
    [a,b]=ismember(Names,inds);
    Names(a)=syms(b(a));
    
    %% Write to File
    data=table(Names,ls{:});
    data.Properties.VariableNames = {'Gene','MeTeOR_Weight','MeTeOR_CN','MeTeOR_AA','BIOGRID','MSigDB','STRINGLIT','EVEX'};
    data = sortrows(data,'MeTeOR_Weight','descend');
    writetable(data,sprintf('%svalidation/Naive%s_%s%s',root,num2str(k),'MeTeOR',strrep(entity,'.','')),'Delimiter','\t');
end

end

