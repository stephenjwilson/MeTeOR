function [NMF_auc,CN_auc,AA_auc] = runGSretro_alt(TestNetName,TestOldNetName,AllNetworks,root,k)
%runGSretro Takes in a present and an old Gold Standard network name as
%well as a test network name (present and past).
%   All strings must be titles that are present in All Networks, root is
%   the central location from which results will be stored. K is the value
%   of factorization
disp(TestNetName);
olddata=AllNetworks.(TestOldNetName);
oldNetwork=olddata{1};
oldNames=olddata{2};

newdata=AllNetworks.(TestNetName);
newNetwork=newdata{1};
newNames=newdata{2};

%disp('NMF');
[result,predNames]=getNMFresult(TestOldNetName,k,AllNetworks,root);

[NMF_auc]=GSretro_FastAlt(oldNames,oldNetwork, newNetwork, newNames, result,predNames);

% run CN and AA
[CN,AA]=naivePrediction(TestOldNetName,AllNetworks);
disp('CN');
[CN_auc]=GSretro_FastAlt(oldNames,oldNetwork, newNetwork, newNames, CN,oldNames);
disp('AA');
[AA_auc]=GSretro_FastAlt(oldNames,oldNetwork, newNetwork, newNames, AA,oldNames);
disp('');
end