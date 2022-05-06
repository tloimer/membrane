% Compare enthalpy differences h(T1,p1) - h(T1, p2) computed by integrating
% dh/dp between p1 and p2, or first integrating along the isenthalpic line,
% finding h(T2,p2) = h(T1, p1) and integrating cp between T2 and T1.

if ~exist('substance.m')
    addpath('~/matlab/program');
end


T1 = 295;
p = 1e5* [...
3.0 3.0 3.0 3.0 3.0 3.0 2.0 2.0 2.0 2.0 1.0 1.0;...
2.5 2.0 1.5 1.0 0.5 0.1 1.5 1.0 0.5 0.1 0.5 0.1];
iso = substance('isobutane');
R142 = substance('R142b');

for i = 1:size(p,2)
    for s = {iso, R142} % s is now a 1x1 cell array
        [T2 hdhdp hcpdT] = compare_enthalpy(s{1}, T1, p(1,i), p(2,i));
        fprintf(['%10s: DTJT = %.4f K, Delta h = %5.0f J, '...
                'difference = %+8.2g\n'],...
                s{1}.name, T1 - T2, hdhdp, (hcpdT-hdhdp)/hdhdp);
    end
end

%%% nested function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T2, hdhdp, hcpdT] = compare_enthalpy(s, T1, p1, p2)

% move out of the loop, for efficiency
%s = substance(name);

% Downstream temperature for isenthalpic process
T2 = s.intjt(T1,p1,p2);

% Consider a T-s diagram,
%     T
%      ^                 p = const
%      | h(T1,p1)      .
%  T1 -|- +----------+ h(T1, p2)
%      |   `       .
%      |     `   .
%  T2 -|-      `
%      |        isenthalpic line
%      +--------------> s
%

hdhdp = -integral(@(p) s.dhdp(T1,p), p2, p1);
hcpdT = integral(@(T) s.cpg(T,p2), T2, T1);
end
