function [sp,baucs,bpr,bpr5,bf1,bf15] = runGSretro(GSNetName,GSoldNetName,TestNetName,TestOldNetName,AllNetworks,root,k)
%runGSretro Takes in a present and an old Gold Standard network name as
%well as a test network name (present and past).
%   All strings must be titles that are present in All Networks, root is
%   the central location from which results will be stored. K is the value
%   of factorization
GSnewdata=AllNetworks.(GSNetName);
GSolddata=AllNetworks.(GSoldNetName);
GSoldNetwork=GSolddata{1};
GSoldNames1=GSolddata{2};
GSoldNames2=GSolddata{2};
GSnewNetwork=GSnewdata{1};
GSnewNames1=GSnewdata{2};
GSnewNames2=GSnewdata{2};

olddata=AllNetworks.(TestOldNetName);
oldNetwork=olddata{1};
oldNames1=olddata{2};
oldNames2=olddata{2};
newdata=AllNetworks.(TestNetName);
newNetwork=newdata{1};
newNames1=newdata{2};
newNames2=newdata{2};

sp={};
c=0;
baucs=zeros(1000,3);
bpr=zeros(1000,3);
bpr5=zeros(1000,3);
bf1=zeros(1000,3);
bf15=zeros(1000,3);
[result,predNames1]=getNMFresult(TestOldNetName,k,AllNetworks,root);

[BAUCs,BPR,BPR5,BF1,BF15]=GSretro(oldNames1,oldNames2,result,predNames1,predNames1, GSoldNetwork, GSoldNames1, GSoldNames2, GSnewNetwork, GSnewNames1,GSnewNames2, sprintf('MeTeOR-%s_NMF%0.1f',TestNetName,k),root);

sp{1}=sprintf('%s',TestOldNetName);
baucs(1:nnz(BAUCs),1)=BAUCs;
bpr(1:nnz(BAUCs),1)=BPR;
bpr5(1:nnz(BAUCs),1)=BPR5;
bf1(1:nnz(BAUCs),1)=BF1;
bf15(1:nnz(BAUCs),1)=BF15;

% run CN and AA
[CN,AA]=naivePrediction(TestOldNetName,AllNetworks);

[BAUCsCN,BPRCN,BPR5CN,BF1CN,BF15CN]=GSretro(oldNames1,oldNames2,CN,predNames1,predNames1, GSoldNetwork, GSoldNames1, GSoldNames2, GSnewNetwork, GSnewNames1,GSnewNames2, sprintf('MeTeOR-%s_CN',TestNetName),root);
[BAUCsAA,BPRAA,BPR5AA,BF1AA,BF15AA]=GSretro(oldNames1,oldNames2,AA,predNames1,predNames1, GSoldNetwork, GSoldNames1, GSoldNames2, GSnewNetwork, GSnewNames1,GSnewNames2, sprintf('MeTeOR-%s_AA',TestNetName),root);


sp{2}=sprintf('%s_CN',TestOldNetName);
baucs(1:nnz(BAUCsCN),2)=BAUCsCN;
bpr(1:nnz(BAUCsCN),2)=BPRCN;
bpr5(1:nnz(BAUCsCN),2)=BPR5CN;
bf1(1:nnz(BAUCsCN),2)=BF1CN;
bf15(1:nnz(BAUCsCN),2)=BF15CN;

sp{3}=sprintf('%s_AA',TestOldNetName);
baucs(1:nnz(BAUCsAA),3)=BAUCsAA;
bpr(1:nnz(BAUCsAA),3)=BPRAA;
bpr5(1:nnz(BAUCsAA),3)=BPR5AA;
bf1(1:nnz(BAUCsAA),3)=BF1AA;
bf15(1:nnz(BAUCsAA),3)=BF15AA;
end