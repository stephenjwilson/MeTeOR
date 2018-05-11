function [ AllNetworks ] = importNetworks(AllNetworks, config)
%importNetworks Crawls the directory of networks and imports them into a
%struct called AllNetworks.
%   This function only imports files with the txt or tsv extensions and
%   that have the word mapped on the end
    netdir=strcat('../',config.general.networkdir);
    D=rdir(strcat(netdir,'**/*.txt'));
    % Defines prefix for each header
    hmap=struct('EntrezID','G.','MeSHID','D.','CID','C.');
    
    for i=1:length(D)
        nm=D(i).name
        origfl=nm;
        nm=strsplit(nm,'/');
        flname=strrep(nm{end},'.txt','');
        if sum(ismember(fields(AllNetworks),flname))==0
            edgelist=parseNetwork(origfl);
            %% Add C,G,D labels for entitiies
            % get header
            fid = fopen(origfl, 'r');
            header = textscan(fid,'%s%s%f%s%[^\n\r]', 1);
            fclose(fid);
            h1=header{1}{1}(2:end); % header label 1
            h2=header{2}{1}; % header label 2
            %% Modify entity names
            edgelist.Entity1=strcat(hmap.(h1),edgelist.Entity1);
            edgelist.Entity2=strcat(hmap.(h2),edgelist.Entity2);
            %% Make symm matrix
            [net,names]=makeSymmMatrix(edgelist,flname,netdir);
            AllNetworks.(flname)={net,names,nm};
        end
    end
    save(strcat(netdir,'AllNetworks'),'AllNetworks','-v7.3')
end

