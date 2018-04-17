%Test MRR

%% Generate Data
testmats= {
    [1 2 0 0; 0 0 3 0; 1 2 2 0]
    [2 1 0 0; 0 0 3 0; 1 2 2 0]
    [0 1 2 0; 0 0 3 0; 1 2 2 0]
    [0 0 1 0; 0 0 3 0; 1 2 2 0]
    [0 1 0 0; 0 0 3 0; 1 2 2 0]};

names={'A','B','C','D'};


%% Test MRR
MRR(testmats, 1,0)