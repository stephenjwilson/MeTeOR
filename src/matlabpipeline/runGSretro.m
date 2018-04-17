function [sp,baucs,bpr,bpr5,bf1,bf15] = runGSretro( GSNetName,GSoldNetName,TestNetName,TestOldNetName,AllNetworks,root,k)
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

ks={k};
%ks={10,30,60,90,120,150,180};

sp={};
c=0;
baucs=zeros(1000,length(ks));
bpr=zeros(1000,length(ks));
bpr5=zeros(1000,length(ks));
bf1=zeros(1000,length(ks));
bf15=zeros(1000,length(ks));
for i=1:length(ks)
    k=ks{i}
    [result,predNames1]=getNMFresult(TestOldNetName,k,AllNetworks,root);
    
    [BAUCs,BPR,BPR5,BF1,BF15]=GSretro(oldNames1,oldNames2,result,predNames1,predNames1, GSoldNetwork, GSoldNames1, GSoldNames2, GSnewNetwork, GSnewNames1,GSnewNames2, sprintf('MeTeOR-%s %0.1f',TestNetName,ks{i}),root);
    [~,l]=size(BAUCs);
    
    for o=1:l
        c=c+1;
        sp{c}=sprintf('%s',TestOldNetName);
        baucs(1:nnz(BAUCs),i)=BAUCs;
        bpr(1:nnz(BAUCs),i)=BPR;
        bpr5(1:nnz(BAUCs),i)=BPR5;
        bf1(1:nnz(BAUCs),i)=BF1;
        bf15(1:nnz(BAUCs),i)=BF15;
    end
end
end