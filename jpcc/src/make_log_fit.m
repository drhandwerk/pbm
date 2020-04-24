function make_log_fit(fit_type,mech,kconst,kvar,ic_in,tend,maxsize,lb,ub,dirID,timeID)

fileID = strcat(dirID,'/',timeID,'_params.log');
f = fopen(fileID,'w');

fprintf(f,'%% Initial values:\n');
fprintf(f,'fit_type = ''%s'';\n',fit_type);
fprintf(f,'mech = ''%s'';\n',mech);

fprintf(f,'kconst = [');
for i = 1:size(kconst)-1
    fprintf(f,'%g; ',kconst(i));
end
fprintf(f,'%g];\n',kconst(end));

fprintf(f,'kvar = [');
for i = 1:size(kvar)-1
    fprintf(f,'%g; ',kvar(i));
end
fprintf(f,'%g];\n',kvar(end));

fprintf(f,'lb = [');
for i = 1:size(lb)-1
    fprintf(f,'%g; ',lb(i));
end
fprintf(f,'%g];\n',lb(end));

fprintf(f,'ub = [');
for i = 1:size(ub)-1
    fprintf(f,'%g; ',ub(i));
end
fprintf(f,'%g];\n',ub(end));

fprintf(f,'ic = %g;\n',ic_in);
fprintf(f,'tend = %g;\n',tend);
fprintf(f,'maxsize = %g;\n',maxsize);
fclose(f);
end

