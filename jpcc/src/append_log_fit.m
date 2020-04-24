function append_log_fit(kargs,dirID,timeID)
fileID = strcat(dirID,'/',timeID,'_fit.log');
f = fopen(fileID,'a');
fprintf(f,'[');
for i = 1:size(kargs)-1
    fprintf(f,'%g; ',kargs(i));
end
fprintf(f,'%g];\n',kargs(end));
fclose(f);
end