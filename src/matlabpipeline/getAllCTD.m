function [names] = getAllCTD(netdir)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
D=rdir(strcat(netdir,'**/CTD_CG*.txt'));
names={};
for i=1:length(D)
    nm=D(i).name;
    nm=strsplit(nm,'/');
    flname=strrep(nm{end},'.txt','');
    names{i}=flname;
end
end

