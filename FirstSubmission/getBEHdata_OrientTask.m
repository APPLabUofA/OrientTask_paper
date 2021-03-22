function [error_deg] = getBEHdata_OrientTask(subj_id, exp)
% Used in the pre-processing code


% Load each subject's data
load([exp.dataLocation '\RawData\BEH\' subj_id '_Orient_v3.mat'])


% Remove practice trials (first 20 trials)
error_deg = data.errorDegrees(1,21:end);
% error_rad = data.errorRads(1,21:end);
% lag_errors = data.lags(1,21:length(data.errorDegrees));
targets = data.target(1,21:length(data.errorDegrees));
% report_deg = data.reportedOrientDegrees(1,21:end);
% report_rad = data.reportedOrientRads(1,21:end);
% report_respdist = data.reportedRespDistance(1,21:end);
% present_deg = data.presentedOrientDegrees(1,21:end);
% present_rads = data.presentedOrientRads(1,21:end);
% diff_deg = data.error_differenceDegrees(1,21:end);
% diff_rads = data.error_differenceRads(1,21:end);

% Remove trials where no target appeared
error_deg = error_deg(targets == 1);
% error_rad = error_rad(targets == 1);
% report_deg = report_deg(targets == 1);
% report_rad = report_rad(targets == 1);


% this corrects a problem with the subj6's recording (there are 2 fewer EEG
% trials than BEH trials)
if strcmpi(subj_id,'subj6')
   error_deg(122:123) = [];
end




















