function [CN,AA,Names1] = naivePrediction(NetName,AllNetworks)
%naivePrediction Summary of this function goes here
%   Detailed explanation goes here
%% Get Data
data=AllNetworks.(NetName);
network=data{1}>0;
Names1=data{2};

try
    load(sprintf('%s_naivepred.mat',NetName), 'CN', 'AA');
catch
    network(logical(eye(length(network))))=0;

    %% Transform the Network by equation
    % Find Common Neighbors
    disp('Setting up naive pred');
    [i,j] = find(tril(ones(size(network)),-1));  % much faster than nchoosek
    CN_values = zeros(length(i),1);
    AA_values = zeros(length(i),1);
    
    degree = sum(network,1);

    disp('Calculating common neighbors');
    % Calculate CN and AA
    parfor ind=1:length(i)
        if mod(ind,100000)==0
            disp(ind);
        end
        common_neighbors = network(i(ind),:)& network(j(ind),:);
        CN_values(ind) = sum(common_neighbors);
        if CN_values(ind)>0
            AA_values(ind) = sum(1 ./ log(degree(common_neighbors)));
        end
    end
    % place into a network
    [m,n] = size(network);
    CN = sparse(i,j,CN_values,m,n);
    AA = sparse(i,j,AA_values,m,n);
    save(sprintf('%s_naivepred.mat',NetName), 'CN', 'AA');
end
end

