function curv=curv()
%CURV       Curvature [1/m].
%
%  Calls COSTHETA, EPSILON, KAPPA.

% parallel plate: kap = eps*b^2/12,
%   pcap = sig/(b/2) = 2*sig/b.
% parallel plates
% pcap = sig(380)/sqrt(3*kappa/epsilon);
% capillary bundle: kap = eps*R^2/6,
%   pcap = 2*sig/R.
% capillary bundle
% pcap = sig(380)*costheta/sqrt(1.5*kappa/epsilon);

% epsilon: Lückengrad - volume fraction of pores

% capillary bundle
%curv = costheta/sqrt(1.5*kappa/epsilon);
% parallel plates
curv = costheta/sqrt(3*kappa/epsilon);
