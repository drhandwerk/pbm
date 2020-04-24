function log_mepbm(mech,kargs,ic_in,tend,maxsize,sol)
% check if directory exists and create if not
sl = './sim_logs';
if (0 == exist(sl,'dir'))
    mkdir(sl);
end
timeID = datestr(now,'yyyy-mm-dd_HH-MM-SS');
dirID = strcat(sl,'/',timeID);
mkdir(dirID);

fileID = strcat(dirID,'/params.log');
f = fopen(fileID,'w');
fprintf(f,'mech = ''%s'';\n',mech);
fprintf(f,'kargs = [');
for i = 1:size(kargs)-1
    fprintf(f,'%g; ',kargs(i));
end
fprintf(f,'%g];\n',kargs(end));
fprintf(f,'ic = %g;\n',ic_in);
fprintf(f,'tend = %g;\n',tend);
fprintf(f,'maxsize = %g;\n',maxsize);
fclose(f);
write_sol(strcat(dirID,'/',timeID,'_sol'),sol);
end

