clear all
clc

load('caseA.mat')
ward_pos_caseA = ward_pos_case;
ward_cov_caseA = ward_cov_case;

load('caseB.mat')
ward_pos_caseB = ward_pos_case;
ward_cov_caseB = ward_cov_case;

load('caseC.mat')
ward_pos_caseC = ward_pos_case;
ward_cov_caseC = ward_cov_case;

load('caseD.mat')
ward_pos_caseD = ward_pos_case;
ward_cov_caseD = ward_cov_case;

load('caseE.mat')
ward_pos_caseE = ward_pos_case;
ward_cov_caseE = ward_cov_case;

load('caseF.mat')
ward_pos_caseF = ward_pos_case;
ward_cov_caseF = ward_cov_case;

load('caseG.mat')
ward_pos_caseG = ward_pos_case;
ward_cov_caseG = ward_cov_case;

save('ward.mat', ...
    "ward_pos_caseA", "ward_cov_caseA", ...
    "ward_pos_caseB", "ward_cov_caseB", ...
    "ward_pos_caseC", "ward_cov_caseC", ...
    "ward_pos_caseD", "ward_cov_caseD", ...
    "ward_pos_caseE", "ward_cov_caseE", ...
    "ward_pos_caseF", "ward_cov_caseF", ...
    "ward_pos_caseG", "ward_cov_caseG")