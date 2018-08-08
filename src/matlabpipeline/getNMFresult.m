function [ result , Names1] = getNMFresult(NetName,k,AllNetworks,root,rep,trimfactor,iter,maxiter)
%getNMFresult Summary of this function goes here
%   [ result , Names1] = getNMFresult(NetName,k,AllNetworks, rep,trimfactor, iter)
%   Takes a network name in the allnetworks struct. Also takes a k for NMF,
%   rep for the number of replicates to do in intialization, trimfactor to
%   eliminate low-confidence edges, and iter to label the current iteration
%   of NMF calculation

if nargin<5
    rep=3;
end
if nargin<6
    trimfactor=2; % gets rid of entities that have less than the specified number of edges
end
if nargin<7
    iter='';
end
if nargin<8
    maxiter=100;
end

data=AllNetworks.(NetName);
Network=data{1};
Names1=data{2};
%Names2=data{2};

%% Start or load NMF
try
    if ~contains(NetName,'MRCOC')
        load(sprintf('%sMeTeOR%sFact%s%s%s',root,num2str(k),num2str(trimfactor),NetName,iter));
    else
        load(sprintf('%s%sFact%s%s%s',root,num2str(k),num2str(trimfactor),NetName,iter));
    end
    %result=sparse(result);
    %save(sprintf('MeTeOR%sFact%s%s_fast',num2str(k),num2str(trimfactor),NetName),'result','Names1','-v7.3');
    disp('Loaded from file');
catch
    disp('File Not Found');
    %rng(10,'multFibonacci');
    %seed = 1;
    %n = 10;
    %% Trim Network
    if ~contains(NetName,'MRCOC')
        Names3=cellfun(@(x) x(2)=='.',Names1);
        Names1=Names1(Names3);
        Network=Network(Names3,Names3);
    end
    
    % Trim by number of associations
    s1=sum(Network>0,1);
    Names1=Names1(s1>trimfactor);
    Network=Network(s1>trimfactor,s1>trimfactor);
    while sum(s1<trimfactor)>0
        s1=sum(Network>0,1);
        Names1=Names1(s1>trimfactor);
        Network=Network(s1>trimfactor,s1>trimfactor);
        s1=sum(Network>0,1);
    end
    %% Start NMF
    opt = statset('MaxIter',5,'Display','final');%'Streams',RandStream.create('mrg32k3a','NumStreams',n,'Seed',seed),'UseSubstreams',true);
    [w,h] = nnmf(Network,k,'algorithm','mult','rep',rep,'options',opt);
    opt = statset('Maxiter',maxiter,'Display','final'); 
    [w,h] = nnmf(Network,k,'w0',w,'h0',h,...
                 'options',opt,...
                 'algorithm','als');
    result=w*h;
    result=(result+result')/2;
    result(result<.1)=0;
    result=sparse(result);
    if ~contains(NetName,'MRCOC')
        save(sprintf('%sMeTeOR%sFact%s%s%s',root,num2str(k),num2str(trimfactor),NetName,iter),'result','Names1','-v7.3');
    else
        save(sprintf('%s%sFact%s%s%s',root,num2str(k),num2str(trimfactor),NetName,iter),'result','Names1','-v7.3');
    end
    clear('w','h');
end
end
