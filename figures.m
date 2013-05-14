% Create the figures for the proceeding contribution
% Copy from figs.m
addpath('program'); % a no-op if 'program' is already added

ispgfplot = true;
%ispgfplot = false;

T1 = 300;
theta = 0;
s = substance('butane');
f = fmodel('homogeneous');
mem = membrane(60e-9,0.44,0.5,'tube',1,8.1,1e-3);

% VAPOR FLOW
p1 = 2.1e5;
p2 = 5e4;

[m fl] = mnum(T1,p1,p2,theta,s,mem,f);

pTplot(fl,'f2pT',ispgfplot);
pTzplots(fl,'f2',ispgfplot);
Tsplot(fl,'f2Ts',ispgfplot);

% PARTIAL CONDENSATION CASE
p1=2.3e5;

[m fl] = mnum(T1,p1,p2,theta,s,mem,f);

pTplot(fl,'f3pT',ispgfplot);
pTzplots(fl,'f3',ispgfplot);
Tsplot(fl,'f3Ts',ispgfplot);

% COMPLETE CONDENSATION

mem = membrane(26e-9,0.44,0.5,'tube',1,8.1,1e-3);

p1 = s.ps(T1);

[m fl] = mnum(T1,p1,p2,theta,s,mem,f);

pTplot(fl,'f4pT',ispgfplot);
pTzplots(fl,'f4',ispgfplot);
Tsplot(fl,'f4Ts',ispgfplot);
