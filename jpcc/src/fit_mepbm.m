function [k,fval,exitflag,output] = fit_mepbm(mech,kconst,kvar,ic_in,tend,maxsize,fit_type,fit_data,lb,ub,options)
if strcmp(fit_type,'precursor') && tend < fit_data(end,1)
    error('Precursor data end time larger than tend.');
end
if nargin < 11
    options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf,'MaxIterations',50,'CompletePoll','On');
end
if nargin < 9
    lb = zeros(length(kvar),1);
    ub = Inf(length(kvar),1);
end

% logging
[dirID,timeID] = fit_log_mkdir();
make_log_fit(fit_type,mech,kconst,kvar,ic_in,tend,maxsize,lb,ub,dirID,timeID);

% make interpolated domain
xq = linspace(0,4.5,1000);
if strcmp(fit_type,'histogram')
        h1 = histfit(fit_data,25,'kernel'); xlim([0 4.5]); ylim([0 Inf]); xlabel('nm'); ylabel('TEM data counts');
        u = h1(2).XData;
        v = h1(2).YData;
        ext_l = min(u):-mean(diff(u)):0;
        ext_l = flip(ext_l);
        ext_l = ext_l(1:end-1);
        ext_r = max(u):mean(diff(u)):5;
        ext_r = ext_r(2:end);
        u_ext = [ext_l u ext_r];
        v_ext = [zeros(1,length(ext_l)) v zeros(1,length(ext_r))];
        v_ext_n = v_ext./max(v_ext);
        F = griddedInterpolant(u_ext,v_ext_n);
        vq = F(xq);
end

switch fit_type
    case 'histogram'
        [k,fval,exitflag,output] = patternsearch(@(kvar) hist_obj(kvar),kvar,[],[],[],[],lb,ub,[],options);
    case 'precursor'
        [k,fval,exitflag,output] = patternsearch(@(kvar) precursor_obj(kvar),kvar,[],[],[],[],lb,ub,[],options);
    otherwise
        error('Unrecognized fit type');
   
end

append_log_params(k,fval,dirID,timeID)

    function d = hist_obj(kvar)
        append_log_fit(kvar,dirID,timeID);
        sol = solve_mepbm(mech,[kconst; kvar],ic_in,tend,maxsize);
        TEMsol = deval(sol,tend);
        TEMsol = TEMsol(3:maxsize);
        TEMsol_ext = [0 0 0 TEMsol' zeros(1,4627-(length(TEMsol)+3))];
        TEMsol4_ext_n = TEMsol_ext./max(TEMsol_ext);
        F = griddedInterpolant(atoms_to_diam(0:4626),TEMsol4_ext_n);
        vqd = F(xq);
        d = trapz(abs(vq-vqd)); % L1 norm!
    end

    function d = precursor_obj(kvar)
        append_log_fit(kvar,dirID);
        sol = solve_mepbm(mech,[kconst; kvar],ic_in,tend,maxsize);
        f = deval(sol,fit_data(:,1));
        d = norm(fit_data(:,2)./1375 - f(1,:)');
    end

end