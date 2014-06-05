if ~exist('substance.m'), addpath('../program'); end

%Setup the flow problem
T1 = 293;
p1 = 1.1e5;
s = substance('isobutane');
f = fmodel('plug');
theta = 0;

pu1 = membrane(10e-9,0.5,24,'tube',2,8.1,20e-6);
pu2 = membrane(100e-9,0.4,24,'tube',2.5,8.1,150e-6);
pu3 = membrane(6e-6,0.3,24,'tube',3,8.1,2e-3);

fl1 = flowsetup(T1,T1,theta,s,pu1,f);
fl2 = flowsetup(T1,T1,theta,s,pu2,f);
fl3 = flowsetup(T1,T1,theta,s,pu3,f);

fl1.knudsen(T1,p1)
fl2.knudsen(T1,p1)
fl3.knudsen(T1,p1)

fprintf('For testing, specific gas const = %g J/kgK, nu = %g m2/s.\n',s.R,s.nug(T1,p1));
