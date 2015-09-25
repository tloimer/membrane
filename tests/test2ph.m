%TEST2PH    Trigger a bug - upstream two-phase mixture

if ~exist('substance.m','file'), addpath('../program'); end
%------------------------------------------------INPUT
%Setup the flow problem
T1 = 293.15;
%s = substance('butane');
s = substance('isobutane');
%s = substance('nitrogen');
%f = fmodel('me2');
f = fmodel('plug');

pu1 = membrane(30e-9,0.8,11.8,'tube',3,8.1,10e-6);
pu2 = membrane(100e-9,0.6,35,'tube',3,8.1,25e-6);
pu3 = membrane(400e-9,0.5,35,'tube',3,8.1,25e-6); %400
pu4 = membrane(800e-9,0.5,30,'tube',3,8.1,25e-6);
pu5 = membrane(2.5e-6,0.4,30,'tube',3,8.1,1e-3);
pu6 = membrane(2.5e-6,0.4,1.13,'tube',3,8.1,1.83e-3);
%----------------------------------------------------------
fprintf([upper(mfilename)...
	':  Test cases with two-phase flow upstream of a membrane\n']);

prmems = {{pu6 pu5 pu4 pu3 pu2 pu1}}; pmr = mstackstruct(0,prmems,f);
mbug2 = {{pu1 pu5},{pu6 pu5 pu4 pu3 pu2 pu1}}; msbug = mstackstruct(0,mbug2,f);
psat1 = s.ps(T1);
deltap = 0.03e5;

m2bug = mnumadiabat(T1,psat1,psat1-deltap,s,msbug);
fprintf('Two membranes, 6 layers + 2 layers: m = %g,\n', m2bug);
mbug = mnumadiabat(T1,psat1,psat1-deltap,s,pmr);
fprintf('One membrane, 6 layers: m = %g,\n', mbug);

% Try to trigger the error, by directly calling asym
m = 0.03276; % m = 0.0327615366317;
T2bug = 293.073; % T1bug = 293.072929543;
p2bug = 298175.; % p2bug = 298175.484228
s2 = downstreamstate(T2bug,p2bug,1,0,s,m);
pmr.T1 = T1; pmr.p1in = psat1; pmr.T2 = T2bug; pmr.p2 = p2bug;
pmr.a = 1; pmr.q2 = 0;
pmr = pmr.writeflowsetups(T1,T2bug,s,pmr);
%[p1asym pmr] = asym(m,s2,pmr,solverstruct('accurate'));
if ~exist('p1asym','var') || isempty(p1asym)
  fprintf('p1asym is empty!\n');
else
  fprintf('p1asym = %.6g\n', p1asym);
end

% Now try to compute the second membrane, with two layers
msbug.T1 = T1; msbug.p1in = psat1; msbug.T2 = T2bug; msbug.p2 = p2bug;
msbug.a = 1; msbug.q2 = 0;
msbug = msbug.writeflowsetups(T1,T2bug,s,msbug);
[p1asym msbug] = asym(m,s2,msbug,solverstruct('accurate'));
if isempty(p1asym)
  fprintf('p1asym is empty!\n');
else
  fprintf('p1asym = %.6g\n', p1asym);
end

clear('-regexp','pu?','T1','s','f','m','prmems','pmr','mbug2','msbug',...
	'psat1','deltap','m2bug','mbug','s2','T2bug','p2bug','p1asym');
fprintf([upper(mfilename) ' finished.\n']);
