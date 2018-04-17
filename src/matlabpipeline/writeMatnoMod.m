function [  ] = writeMatnoMod( flname,titles,mat )
%writeMatnoMod Helper function to write data to file
%   Writes data to csv files

fid = fopen(strcat(flname,'.csv'), 'w') ;
fprintf(fid, '%s,', titles{:}) ;
fprintf(fid, '\n');
fclose(fid);
dlmwrite(strcat(flname,'.csv'),mat,'-append');

fid = fopen(strcat(flname,'summary.csv'), 'w') ;
fprintf(fid, ',%s,,,\n','MeTeOR') ;
%fprintf(fid, '%s,', titles{:}) ;
fprintf(fid, ',mean,SD,N\n');
[~,c]=size(mat);
%me=mean(mat);
%sd=std(mat);



f='%f,%f,%f\n';
for o=1:c
    fprintf(fid,'%s,',titles{o});
    dat=mat(:,o);
    dat=dat(dat>0);
    me=mean(dat);
    sd=std(dat);
    N=sum(dat>0);
    fprintf(fid, f, me,sd,N);
    sprintf(f, me,sd,N)
   
end


end

