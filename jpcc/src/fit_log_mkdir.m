function [dirID,timeID] = fit_log_mkdir()
% check if directory exists and create if not
fitDir = './fit_logs';
if (0 == exist(fitDir,'dir'))
    mkdir(fitDir);
end
timeID = datestr(now,'yyyy-mm-dd_HH-MM-SS');
dirID = strcat(fitDir,'/',timeID);
mkdir(dirID);

end