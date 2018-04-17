%% Experimenting with Network Characterization

network=AllNetworks.MeTeORallMeTeORallGeneGene{1};
Names1=AllNetworks.MeTeORallMeTeORallGeneGene{2};

% Number of non-zero
disp('The number of non-zero elements')
nnz(network)

% Sumations
s1=sum(network,1);
s2=sum(network,2);

disp('The sums of the dimensions are equal')
sum(s1==s2')
disp('The sums of the dimension are 1')
sum(s1==1)
disp('The sums of the dimension are 0')
sum(s2==0)

%% Trim the network
disp('The trimmed network')
Names2=cellfun(@(x) x(2)=='.',Names1);
Names3=Names1(Names2);
newNetwork=network(Names2,Names2);

% Number of non-zero
disp('The number of non-zero elements')
nnz(newNetwork)

% Sumations

s1=sum(newNetwork,1);
s2=sum(newNetwork,2);

disp('The sums of the dimensions are equal')
sum(s1==s2')
disp('The sums of the dimension are 1')
sum(s1==1)
disp('The sums of the dimension are 0')
sum(s2==0)

%% Improved Trimmed network
disp('The Improved Trimmed network')
Names2=cellfun(@(x) x(2)=='.',Names1);
Names4=Names1(Names2);
newNetwork=network(Names2,Names2);

% Trim by number of associations
s1=sum(newNetwork,1);
Names4=Names4(s1>1);
newNetwork2=newNetwork(s1>1,s1>1);

while sum(s1==1)>0
    s1=sum(newNetwork2,1);
    Names4=Names4(s1>1);
    newNetwork2=newNetwork2(s1>1,s1>1);
end

% Number of non-zero
disp('The number of non-zero elements')
nnz(newNetwork2)

% Sumations
s1=sum(newNetwork2,1);
s2=sum(newNetwork2,2);

disp('The sums of the dimensions are equal')
sum(s1==s2')
disp('The sums of the dimension are 1')
sum(s1==1)
disp('The sums of the dimension are 0')
sum(s2==0)

