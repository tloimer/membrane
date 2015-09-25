function testcond()
%TESTCOND   A simple test of the conductivity models.

if ~exist('substance.m','file'), addpath('../program'); end
%------------------------------------------------INPUT
fprintf([upper(mfilename)...
	'...run all conductivity models on an arbitrary test case.\n']);
%Setup the flow problem
T1 = 293.15;
%s = substance('butane');
s = substance('isobutane');
%s = substance('nitrogen');
psat1 = s.ps(T1);
deltap = 0.03e5;

pu1 = membrane(30e-9,0.8,11.8,'tube',3,8.1,10e-6);
pu2 = membrane(100e-9,0.6,35,'tube',3,8.1,25e-6);
pu3 = membrane(400e-9,0.5,35,'tube',3,8.1,25e-6); %400
pu4 = membrane(800e-9,0.5,30,'tube',3,8.1,25e-6);
pu5 = membrane(2.5e-6,0.4,30,'tube',3,8.1,1e-3);
pu6 = membrane(2.5e-6,0.4,1.13,'tube',3,8.1,1.83e-3);
prmems = {{pu6 pu5 pu4 pu3 pu2 pu1}};

cn = {'series','parallel','me1','me2','EMT'};
for i = 1:size(cn,2)
  f = fmodel(cn{i});
  pmr = mstackstruct(0,prmems,f);
  [m,pmr] = mnumadiabat(T1,psat1,psat1-deltap,s,pmr);
  fprintf('  %-9s m = %.5g kg/m2s, p_residual = %.3f Pa\n',...
    [cn{i} ':'], m, pmr.p1sol-pmr.p1in);
end
fprintf([upper(mfilename) ' finished.\n']);
