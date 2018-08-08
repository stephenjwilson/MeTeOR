function [CN,AA,Names1] = naivePrediction(NetName,AllNetworks)
%naivePrediction Summary of this function goes here
%   Detailed explanation goes here
%% Get Data
data=AllNetworks.(NetName);
network=double(data{1}>0);
Names1=data{2};

try
    load(sprintf('%s_naivepred.mat',NetName), 'CN', 'AA');
catch
    network(logical(eye(length(network))))=0;

    %% Transform the Network by equation
    % Find Common Neighbors
    disp('Setting up naive pred');
    CN = network*network;
    mask = eye(length(CN))>0;
    degree = CN(mask);
    CN(mask) = 0;
    AA = diagsum(1./log(degree.*outer(network,network,0)),2,4);
    AA(mask) = 0;
    
    save(sprintf('%s_naivepred.mat',NetName), 'CN', 'AA');
end
end

