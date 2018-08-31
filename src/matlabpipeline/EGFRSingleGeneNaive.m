NetName='MeTeORgenegene';
data=AllNetworks.(NetName);
network=data{1}>0;
Names1=data{2};
entity='G.1956';
ind=find(strcmp(entity,Names1));
network(logical(eye(length(network))))=0;

%% Transform the Network by equation
% Find Common Neighbors
disp('Setting up naive pred');

degree = sum(network,1);

disp('Calculating common neighbors');
% Calculate CN and AA

CN_values = zeros(length(network),1);
AA_values = zeros(length(network),1);
parfor i=1:length(network)
    if mod(i,10000)==0
        disp(i);
    end
    common_neighbors = network(i,:)& network(ind,:);
    CN_values(i) = sum(common_neighbors);
    if CN_values(i)>0
        AA_values(i) = sum(1 ./ log(degree(common_neighbors)));
    end
end

save(sprintf('%s_EGFR.mat',NetName), 'CN_values', 'AA_values');

meteor = AllNetworks.MeTeORgenegene{1}(ind,:);

forwardRank=tiedrank(meteor')'; % forwardRank
meteor_rank=max(forwardRank)-forwardRank+1; % Method Rank

forwardRank=tiedrank(CN_values'); % forwardRank
rank_method=max(forwardRank)-forwardRank+1; % Method Rank
CN_diff=meteor_rank-rank_method; % Relative diff rank

forwardRank=tiedrank(AA_values'); % forwardRank
rank_method=max(forwardRank)-forwardRank+1; % Method Rank
AA_diff=meteor_rank-rank_method; % Relative diff rank

data = table(tmpNames,meteor', AA_values,CN_values, AA_diff', CN_diff');
data.Properties.VariableNames = {'Gene', 'MeTeOR', 'AA', 'CN', 'AADiff','CNDiff'};
data = sortrows(data,'MeTeOR','descend');
writetable(data, '/lab/cedar/home/wilson/tmp_AA.txt','Delimiter','\t');