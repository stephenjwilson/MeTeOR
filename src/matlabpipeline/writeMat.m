function [  ] = writeMat( flname,titles,mat )
%writeMat Helper function to write data to file
%   Writes data to csv files
MOD=2;
fid = fopen(strcat(flname,'.csv'), 'w') ;
fprintf(fid, '%s,', titles{:}) ;
fprintf(fid, '\n');
fclose(fid);
dlmwrite(strcat(flname,'.csv'),mat,'-append');

fid = fopen(strcat(flname,'summary.csv'), 'w') ;
fprintf(fid, ',%s,,,%s,,\n','MeTeOR','Literature Derived') ;
%fprintf(fid, '%s,', titles{:}) ;
fprintf(fid, ',mean,SD,N,mean,SD,N\n');
[~,c]=size(mat);
%me=mean(mat);
%sd=std(mat);
for i=1:MOD:c
	fprintf(fid,'%s,',titles{i});
	f='%f,%f,%f,';
	for o=0:MOD-1
		dat=mat(:,i+o);
		dat=dat(dat>0);
		me=mean(dat);
		sd=std(dat);
		N=sum(dat>0);
		fprintf(fid, f, me,sd,N);
		sprintf(f, me,sd,N)
		f='%f,%f,%f\n';
	end

end

end

