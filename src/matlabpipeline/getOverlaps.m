function [] = getOverlaps(AllNetworks,YamlStruct,test)
%getOverlaps Determines overlaps of the networks of interest
%   Detailed explanation goes here

if nargin<3
    test='MeTeOR';
end

%% Set-Up

% load AllNetworks
if strcmp(test,'MeTeOR')
    titles={{'BIOGRIDCocrystalStructure','BIOGRIDLow','MeTeORgenegene'},{'BIOGRIDchem','DGIdb','MeTeORgenechemical'},{'CTD_human_fromDisgenet','UNIPROT_fromDisgenet','MeTeORdiseasegene'}};
elseif strcmp(test,'Lit')
    titles={{'BIOGRIDCocrystalStructure','BIOGRIDLow','STRING10_textmining'},{'BIOGRIDchem','DGIdb','STITCH5_textmining'},{'CTD_human_fromDisgenet','UNIPROT_fromDisgenet','BeFree'}};
else
    msgID = 'getOverlaps:BadTestCase';
    msg = 'Unable to identify test case. Use MeTeOR or Lit.';
    baseException = MException(msgID,msg);
    throw(baseException);
end


resultdir=strcat('../',YamlStruct.general.resultsdir,'validation/Fig2/overlaps/');
if ~exist(resultdir, 'dir')
    mkdir(resultdir);
end
%% Compare Networks
data={};
for c=1:length(titles)
    GStitles=titles{c};
    comb=combnk(1:length(GStitles),2);
    sub={};
    for o=1:length(comb)
        ind1=comb(o,1);
        ind2=comb(o,2);
        ti=GStitles(ind1);
        tj=GStitles(ind2);

        netname=GStitles{ind1};
       
        GS=AllNetworks.(netname){1};    

        GSnames=AllNetworks.(netname){2};
        
        [w,h]=size(GSnames);
        if h==1
            GSnames1=GSnames;
            GSnames2=GSnames;
        else
            GSnames1=GSnames{1};
            GSnames2=GSnames{2};
        end

        [GS,GSnames1,GSnames2]=eliminateDuplicateNames(GS,GSnames1,GSnames2);

        if length(GSnames1)==length(GSnames2)
           GS=tril(GS);
        end

        netname=GStitles{ind2};

        GS2=AllNetworks.(netname){1};    
        GSnames=AllNetworks.(netname){2};

        [w,h]=size(GSnames);
        if h==1
            GS2names1=GSnames;
            GS2names2=GSnames;
        else
            GS2names1=GSnames{1};
            GS2names2=GSnames{2};
        end

        [GS2,GS2names1,GS2names2]=eliminateDuplicateNames(GS2,GS2names1,GS2names2);

        [GS1_to_common_mapping1]=ismember(GSnames1,GS2names1);
        commonNames1=GSnames1(GS1_to_common_mapping1);
        [GS2_to_common_mapping1]=ismember(GS2names1,commonNames1);

        [GS1_to_common_mapping2]=ismember(GSnames2,GS2names2);
        commonNames2=GSnames2(GS1_to_common_mapping2);
        [GS2_to_common_mapping2]=ismember(GS2names2,commonNames2);

        if length(GS2names1)==length(GS2names2)
           GS2=tril(GS2);
        end
        
        Alli=nnz(GS);
        Allj=nnz(GS2);
        
        ci=nnz(GS(GS1_to_common_mapping1,GS1_to_common_mapping2)>0);
        cj=nnz(GS2(GS2_to_common_mapping1,GS2_to_common_mapping2)>0);
        sub{ismember(GStitles,ti),1}=Alli;
        sub{ismember(GStitles,tj),1}=Allj;
        sub{ismember(GStitles,ti),2}=ti;
        sub{ismember(GStitles,tj),2}=tj;
        
        overlap=nnz(and(GS(GS1_to_common_mapping1,GS1_to_common_mapping2)>0,GS2(GS2_to_common_mapping1,GS2_to_common_mapping2)>0));
        sub{o+3,1}=overlap;
        sub{o+3,2}=ti;
        sub{o+3,3}=tj;
    end
    
    
    
    %% Three way overlap
    netname=GStitles{1};

    GS=AllNetworks.(netname){1};
    GSnames=AllNetworks.(netname){2};
    
    [w,h]=size(GSnames);
    if h==1
        GSnames1=GSnames;
        GSnames2=GSnames;
    else
        GSnames1=GSnames{1};
        GSnames2=GSnames{2};
    end

    [GS,GSnames1,GSnames2]=eliminateDuplicateNames(GS,GSnames1,GSnames2);
    
    if length(GSnames1)==length(GSnames2)
        GS=tril(GS);
    end

    netname=GStitles{2};

    GS2=AllNetworks.(netname){1};
    GSnames=AllNetworks.(netname){2};
    
    [w,h]=size(GSnames);
    if h==1
        GS2names1=GSnames;
        GS2names2=GSnames;
    else
        GS2names1=GSnames{1};
        GS2names2=GSnames{2};
    end
    
    [GS2,GS2names1,GS2names2]=eliminateDuplicateNames(GS2,GS2names1,GS2names2);
    
    if length(GS2names1)==length(GS2names2)
        GS2=tril(GS2);
    end

    netname=GStitles{3};
    GS3=AllNetworks.(netname){1};
    GSnames=AllNetworks.(netname){2};
    
    [w,h]=size(GSnames);
    if h==1
        GS3names1=GSnames;
        GS3names2=GSnames;
    else
        GS3names1=GSnames{1};
        GS3names2=GSnames{2};
    end
    
    [GS3,GS3names1,GS3names2]=eliminateDuplicateNames(GS3,GS3names1,GS3names2);
    
    
    if length(GS3names1)==length(GS3names2)
        GS3=tril(GS3);
    end

    %% names 1 mapping to common
    [GS1_to_common_mapping1]=ismember(GSnames1,GS2names1);
    commonNames=GSnames1(GS1_to_common_mapping1);
    [GS3_to_common_mapping1]=ismember(GS3names1,commonNames);
    
    commonNames1=GS3names1(GS3_to_common_mapping1);
    [GS1_to_common_mapping1]=ismember(GSnames1,commonNames1);
    [GS2_to_common_mapping1]=ismember(GS2names1,commonNames1);
    %% names 2
    [GS1_to_common_mapping2]=ismember(GSnames2,GS2names2);
    commonNames=GSnames2(GS1_to_common_mapping2);
    [GS3_to_common_mapping2]=ismember(GS3names2,commonNames);
    
    commonNames2=GS3names2(GS3_to_common_mapping2);
    [GS1_to_common_mapping2]=ismember(GSnames2,commonNames2);
    [GS2_to_common_mapping2]=ismember(GS2names2,commonNames2);

    overlap=nnz(and(and(GS(GS1_to_common_mapping1,GS1_to_common_mapping2)>0,GS2(GS2_to_common_mapping1,GS2_to_common_mapping2)>0),GS3(GS3_to_common_mapping1,GS3_to_common_mapping2)>0));
    sub{7,1}=overlap;
    sub{7,2}='Overlap';
    data{c}=sub;
end
fh=findall(0,'type','figure');
for i=1:length(fh)
clf(fh(i));
end

%% Make Venn Diagrams and save a summary table
labels={};
vectors=[];
new=[];
for i=1:length(data)

    vect=[data{i}{:,1}];
    vect(1)=vect(1)-vect(4)-vect(5)+vect(7);
    vect(2)=vect(2)-vect(4)-vect(6)+vect(7);
    vect(3)=vect(3)-vect(5)-vect(6)+vect(7);
    vect(4)=vect(4)-vect(7);
    vect(5)=vect(5)-vect(7);
    vect(6)=vect(6)-vect(7);
    
    st=zeros(7,3);
    st(1:3,1)=vect(1:3);
    st(4:6,2)=vect(4:6);
    st(7,3)=vect(7);

    vectors=[vectors vect];
    new=[new;st];
    labels=[labels cellfun(@(x,y) sprintf('%s %s',char(x),char(y)),{data{i}{:,2}},{data{i}{:,3}}, 'UniformOutput',0)];
end
%overlaps=table(labels',vectors');
overlaps=table(labels',new);

save(sprintf('%s%soverlaps.mat',resultdir,test),'overlaps','-v7.3');
writetable(overlaps,sprintf('%s%soverlaps.txt',resultdir,test));

end

