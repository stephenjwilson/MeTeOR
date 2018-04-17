function [mat] = thresholdMat(mat,thresh)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [r,c,v]=find(mat);
    r=r(v<thresh);
    c=c(v<thresh);
    mat(r,c)=0;
end

