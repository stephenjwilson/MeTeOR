function [ reduced_vector ] = reduceVector( vector, size )
%REDUCEVECTOR sample from vector in a deterministic way, to reduce the
%amount of data needed to keep saved
%   gets size number of points from vector, evenly spaced

    if(length(vector) < size)
        reduced_vector=vector;
        while(length(reduced_vector) < size)
            reduced_vector(end+1)=reduced_vector(end);
        end
        
        return;
    end
    
    reduced_vector=zeros(size,1);
    reduced_vector(1)=vector(1);
    count = 2;
    for fraction=[1/(size-1):1/(size-1):1];
        reduced_vector(count)=vector(floor(length(vector)*fraction));
        count = count + 1;
    end
end

