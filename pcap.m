function pcap=pcap(T)
%PCAP(T)    Capillary pressure [Pa].
%
%  Calls SIG.

% parallel plate: kap = eps*b^2/12,
%   pcap = sig/(b/2) = 2*sig/b.
% parallel plates
% pcap = sig(380)/sqrt(3*kappa/epsilon);
% capillary bundle: kap = eps*R^2/6,
%   pcap = 2*sig/R.
% capillary bundle
% pcap = sig(380)/sqrt(1.5*kappa/epsilon);

%pcap=0;return
% epsilon: Lückengrad - volume fraction of pores
epsilon=0.6;

% parallel plates
pcap = sig(380)/sqrt(3*kappa/epsilon);
