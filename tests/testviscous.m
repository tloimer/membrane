function testviscous()
%TESTVISCOUS Test mviscous().
%  Uncomment the line "fprintf('Mass fluex guessed m...)" bevor running this
%  test.

if ~exist('substance.m'), addpath('../program'); end

%Setup the flow problem
T1 = 293.15;
%s = substance('butane');
s = substance('isobutane');
%s = substance('nitrogen');
f = fmodel('plug');

pu1 = membrane(10e-9,0.6,36,'tube',3,8.1,20e-6);
pu2 = membrane(100e-9,0.6,36,'tube',3,8.1,150e-6);
pu3 = membrane(6e-6,0.6,36,'tube',3,8.1,2e-3);

psmems = {{pu1 pu2 pu3}}; pms = mstackstruct(0,psmems,f); pmsorig = pms;
prmems = {{pu3 pu2 pu1}}; pmr = mstackstruct(0,prmems,f); pmrorig = pmr;

psat1 = s.ps(T1);

% Plot with p1 - p2 = 0.1 bar.

%poben = [2.1:0.01:psat1/1e5 psat1/1e5]*1e5;
poben = [1.1:0.2:psat1/1e5 psat1/1e5]*1e5;
mf = poben; mr = poben;

%% mviscous does not need to call ms.writeflowsetup, but it needs the substance
%pms.substance = s;
%pmr.substance = s;

deltap = 0.1e5;
for i = 1:length(poben)
  p2 = poben(i) - deltap;
  mf(i) = mviscous(T1,poben(i),p2,s,pms);
  mr(i) = mviscous(T1,poben(i),p2,s,pmr);
  mguessf = pms.mfluxviscous(T1,poben(i),p2,s,pms);
  mguessr = pmr.mfluxviscous(T1,poben(i),p2,s,pmr);
  fprintf('Mass flux guessed mguess = %8.3g, calculated m = %8.3g, (mguess - m)/m = %.3g%%\n',...
	  mguessf, mf(i), (mguessf - mf(i))*100 / mf(i) );
  fprintf('Mass flux guessed mguess = %8.3g, calculated m = %8.3g, (mguess - m)/m = %.3g%%\n',...
	  mguessr, mr(i), (mguessr - mr(i))*100 / mr(i) );
end

% Plot with p1 - p2 = 0.5 bar.

%poben = [2.1:0.01:psat1/1e5 psat1/1e5]*1e5;
poben = [1.1:0.2:psat1/1e5 psat1/1e5]*1e5;
mf = poben; mr = poben;

deltap = 0.5e5;
for i = 1:length(poben)
  p2 = poben(i) - deltap;
  mf(i) = mviscous(T1,poben(i),p2,s,pms);
  mr(i) = mviscous(T1,poben(i),p2,s,pmr);
  mguessf = pms.mfluxviscous(T1,poben(i),p2,s,pms);
  mguessr = pmr.mfluxviscous(T1,poben(i),p2,s,pmr);
  fprintf('Mass flux guessed mguess = %8.3g, calculated m = %8.3g, (mguess - m)/m = %.3g%%\n',...
	  mguessf, mf(i), (mguessf - mf(i))*100 / mf(i) );
  fprintf('Mass flux guessed mguess = %8.3g, calculated m = %8.3g, (mguess - m)/m = %.3g%%\n',...
	  mguessr, mr(i), (mguessr - mr(i))*100 / mr(i) );
end

%poben = [2.1:0.01:psat1/1e5 psat1/1e5]*1e5;
poben = [1.1:0.2:psat1/1e5 psat1/1e5]*1e5;
mf = poben; mr = poben;

deltap = 1e5;
for i = 1:length(poben)
  p2 = poben(i) - deltap;
  mf(i) = mviscous(T1,poben(i),p2,s,pms);
  mr(i) = mviscous(T1,poben(i),p2,s,pmr);
  mguessf = pms.mfluxviscous(T1,poben(i),p2,s,pms);
  mguessr = pmr.mfluxviscous(T1,poben(i),p2,s,pmr);
  fprintf('Mass flux guessed mguess = %8.3g, calculated m = %8.3g, (mguess - m)/m = %.3g%%\n',...
	  mguessf, mf(i), (mguessf - mf(i))*100 / mf(i) );
  fprintf('Mass flux guessed mguess = %8.3g, calculated m = %8.3g, (mguess - m)/m = %.3g%%\n',...
	  mguessr, mr(i), (mguessr - mr(i))*100 / mr(i) );
end
