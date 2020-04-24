function append_log_params(kargs,bfv,dirID,timeID)
fileID = strcat(dirID,'/',timeID,'_params.log');
f = fopen(fileID,'a');
fprintf(f,'\n%% Results from fit:\n');
fprintf(f,'kargs = [');
for i = 1:size(kargs)-1
    fprintf(f,'%g; ',kargs(i));
end
fprintf(f,'%g];\n',kargs(end));
fprintf(f,'bvf = %g\n',bfv);
fclose(f);
end