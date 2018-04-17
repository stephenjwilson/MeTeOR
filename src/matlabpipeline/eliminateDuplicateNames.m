function [ trimmat,trim1,trim2] = eliminateDuplicateNames( mat,names1,names2, upper_ )
%eliminateDuplicateNames Eliminates all duplicate names as well as blank
%labels
%   Eliminates names if not unique. Non-unique names have their vectors
%   combined through the maximum
    if nargin<4
        upper_=0;
    end
    
    if upper_
        names1=upper(names1);
        names2=upper(names2); 
    end
    
    [trim1,i_names1]=unique(names1);
    [trim2,i_names2]=unique(names2);
    
    %% Names 1
    % Identify duplicates for names2
    duplicate_ind1 = setdiff(1:length(names1), i_names1);
    duplicate_value1 = unique(names1(duplicate_ind1));
    
    % Propogate information across duplicates for names1
    for i=1:length(duplicate_value1)
            ind=strcmp(names1,duplicate_value1{i});
            [~,crossind]=find(mat(ind,:));
            newind=zeros(length(names2),1);
            newind(crossind)=1;
            vals=max(mat(ind,newind>0),[],1); % Max Function is here           
            mat(ind,newind>0)=repmat(vals,length(find(ind)),1);
    end    
    
    %% Names 2
    % Identify duplicates for names2
    duplicate_ind2 = setdiff(1:length(names2), i_names2);
    duplicate_value2 = unique(names2(duplicate_ind2));
    
    % Propogate information across duplicates for names2
    for i=1:length(duplicate_value2)
            ind=strcmp(names2,duplicate_value2{i});
            [crossind,~]=find(mat(:,ind));
            newind=zeros(length(names1),1);
            newind(crossind)=1;
            vals=max(mat(newind>0,ind),[],2); % Max Function is here
            mat(newind>0,ind)=repmat(vals,1,length(find(ind)));
    end    

    %% Prepare Output
    % Trim mat
    trimmat=mat(i_names1,i_names2);
    
    % Eliminate empty labels
    tind1=strcmp('',trim1);
    tind2=strcmp('',trim2);
    trim1(tind1)=[];
    trim2(tind2)=[];
    trimmat=trimmat(~tind1,~tind2);
end

