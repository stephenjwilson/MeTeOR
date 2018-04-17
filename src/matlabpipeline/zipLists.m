function [ out ] = zipLists( inputs)
%zipLists zips lists together
%   Takes cell arrays and zips them together
out=[];
for i=1:length(inputs{1}(1,:))
    for o=1:length(inputs)
        out=[out , inputs{o}(:,i)];
    end
end

end

