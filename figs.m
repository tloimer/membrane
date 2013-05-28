% Create p-T, T-z, p-z and T-s plots
addpath('program'); % a no-op if 'program' is already added

ispgfplot = true;
%ispgfplot = false;

if ~exist('fl')
T1 = 300;  %25 + 273.15;
%theta = 120; % contact angle, degree
% theta = 120;	% NON-WETTING CASE
theta = 0;	% PARTIAL CONDENSATION CASE
% pressure conditions, below

%%%	SUBSTANCE		%%%
s = substance('butane');
% Makes accessible substance properties, e.g., s.ps(T), s.kl(T), ...

%%%	MEMBRANE PROPERTIES 	%%%
% mem = membrane(poredia,eps,km,model,tau,beta,L);
% poredia - pore diameter, eps - void fraction, km - thermal conductivity,
% model - pore topology (porousround, channel, tube), tau - tortuosity,
% beta - molecular flow correction factor, L - tube length
% k [w/mK]: Teflon, 0.26; Glas, 0.5 ( bis -1.38, silica);
% PARTIAL CONDENSATION CASE, VAPOR FLOW
%mem = membrane(60e-9,0.44,0.5,'tube',1,8.1,1e-3);
% NON-WETTING CASE
%mem = membrane(1200e-9,.44,100,'tube',1,8.1,1e-3);
% COMPLETE CONDENSATION
mem = membrane(26e-9,0.44,0.5,'tube',1,8.1,1e-3);
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
p12 = 1.5e5;
%p12 = p1-p2;
%           p1 - pK1 + n*p12
% With P1 = -----------------, we define P10 = pK1 - n*p12.
%           ps1 - pK1 + n*p12
% in mnum>flcalcvars, calc.n = s.jt(T1,p1)*dpK1, not s.jt(T1,ps1).
P10 = pK1 - s.jt(T1,ps1)*dpK1*p12;
% Now all pressures are known and a suitable p1 can be given.
%p1 = ps1 -20*(pK1-ps1);
%p1 = pK1 ;%-(ps1-pK1); %P10 + 0.5*(ps1-pK1);
%p1 = P10;% - 0.1*(ps1-P10);
% VAPOR FLOW
p1 = 2.1e5; p2 = 5e4; p12 = p1 - p2;
%p1 = P10 - 0.1*(ps1-P10)
% PARTIAL CONDENSATION CASE
%p1=2.3e5; p2 = 5e4; p12 = p1 - p2;
%p1 = P10 + 0.25*(ps1-P10); p2 = 4e4;
%          0.2, 0.3 sowie p12 = 1.5e5 ist auch OK; p1 = 1.906e5;
%p1 = P10 + 0.25*(ps1-P10); p2 = p1 - p12;
% NON-WETTING CASE
%p1 = 2.4e5; p2 = 0.4e5; p12 = p1 - p2;
% COMPLETE CONDENSATION
p1 = ps1; p2 = 5e4; p12 = p1 - p2;

%%%	SANITY CHECK	%%%
if p1 - p12 < 0
  error('Pressure difference p12 too large, p2 < 0. p1 = %.0f Pa', p1);
end
%p2 = p1 - p12;

%%%	CALCULATE	%%%
%m = mlinearnew(p1,p2,T1,theta,s,mem,f);
[m fl] = mnum(T1,p1,p2,theta,s,mem,f); %,T1,m);

end

%%%	PLOT		%%%
pTplot(fl,'pT',ispgfplot);
pTzplots(fl,'',ispgfplot);
Tsplot(fl,'Ts',ispgfplot);
