% Create p-T, T-z, p-z and T-s plots
addpath('program'); % a no-op if 'program' is already added

if ~exist('fl')
T1 = 25 + 273.15;
p12 = 1.0e5;
ispgfplot = false;
theta = 0; % contact angle, degree
% pressure conditions, below

%%%	SUBSTANCE		%%%
s = substance('butane');
% Makes accessible substance properties, e.g., s.ps(T), s.kl(T), ...

%%%	MEMBRANE PROPERTIES 	%%%
% mem = membrane(poredia,eps,km,model,tau,beta,L);
% poredia - pore diameter, eps - void fraction, km - thermal conductivity,
% model - pore topology (porousround, channel, tube), tau - tortuosity,
% beta - molecular flow correction factor, L - tube length
mem = membrane(10e-9,1,s.kl(T1),'tube',1,8.1,1e-3);
% eps = 1 ... no effect of km; in addition, km = s.kl(T1)
% tau = 1, beta = 8.1 ... ideal values

% Two-phase flow model. Only 'homogeneous' is available.
f = fmodel('homogeneous');

%%%	BOUNDARY CONDITIONS	%%%
% Determine a T2 by expanding the substance to nearly complete vacuum. This is
% necessary because p1 is not yet known. Complete vacuum is  not possible. 
ps1 = s.ps(T1);
T2 = s.intjt(T1,ps1,10);
flsetup = flowsetup(T2,T1,theta,s,mem,f);
[dpK1 pK1] = flsetup.dpkdT(T1);
%           p1 - pK1 + n*p12
% With P1 = -----------------, we define P10 = pK1 - n*p12.
%           ps1 - pK1 + n*p12
% in mnum>flcalcvars, calc.n = s.jt(T1,p1)*dpK1, not s.jt(T1,ps1).
P10 = pK1 - s.jt(T1,ps1)*dpK1*p12;
% Now all pressures are known and a suitable p1 can be given.
p1 = pK1;

%%%	SANITY CHECK	%%%
if p1 - p12 < 0
  error('Pressure difference p12 too large, p2 < 0. p1 = %.0f Pa', p1);
end
p2 = p1 - p12;

%%%	CALCULATE	%%%
m = mlinearnew(p1,p2,T1,theta,s,mem,f);
[m fl] = mnum(T1,p1,p1-p12,theta,s,mem,f); %,T1,m);

end

%%%	PLOT		%%%
pTplot('pT',fl,ispgfplot);
% flatten z, once
pTzplots('z',fl,ispgfplot);
Tsplot('Ts',fl,ispgfplot);
%Tsnew(fl);
