function testgaseous()
%TESTGASEOUS Test mgaseous().

% Probably derived from ../results/131014/figures.m

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

deltap = 0.1e5;
%poben = [2.1:0.01:psat1/1e5 psat1/1e5]*1e5;
poben = [1.1:0.2:psat1/1e5 psat1/1e5]*1e5;
computeandplot(deltap,poben,T1,s,pms,pmr);

deltap = 0.5e5;
%poben = [2.1:0.01:psat1/1e5 psat1/1e5]*1e5;
poben = [1.1:0.2:psat1/1e5 psat1/1e5]*1e5;
computeandplot(deltap,poben,T1,s,pms,pmr);

deltap = 1e5;
%poben = [2.1:0.01:psat1/1e5 psat1/1e5]*1e5;
poben = [1.1:0.2:psat1/1e5 psat1/1e5]*1e5;
computeandplot(deltap,poben,T1,s,pms,pmr);


%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

function [msv,msk,vis,kn,gas] = compute(deltap,poben,T1,s,ms);
% COMPUTE   Compute mass fluxes for viscous, knudsen and gaseous flow.

% Allocate memory.
len = length(poben);
msv(len) = 0;
msk(len) = 0;
vis(len) = 0;
kn(len) = 0;
gas(len) = 0;

for i = 1:len
  p2 = poben(i) - deltap;
  gas(i) = mgaseous(T1,poben(i),p2,s,ms);
  vis(i) = mgaseous(T1,poben(i),p2,s,ms,'viscous');
  kn(i) = mgaseous(T1,poben(i),p2,s,ms,'knudsen');
  msv(i) = ms.mfluxviscous(T1,poben(i),p2,s,ms);
  msk(i) = ms.mfluxknudsen(T1,poben(i),p2,s,ms);
end


function mfluxplot(tstr,poben,msv,msk,vis,kn,gas);
%MFLUXPLOT  Plot guessed and computed mass fluxes.

poben = poben*1e-5;
h = figure('Name',tstr(tstr ~= '\'));
plot(poben,msk,'rd',poben,kn,'r+',poben,msv,'ks',poben,vis,'kx',...
     poben,gas,'ko',poben,kn+vis,'k.');
legend('knudsen, mstack.','knudsen, calc.', ...
	'viscous, mstack.','viscous, calc.', 'full gaseous', 'kn. + visc., c.');
legend('boxoff');
xlabel('p_1 [bar]');
ylabel('massflux [kg/m2s]');
title(tstr);


function mratioplots(pstr,poben,gas,msv,vis,msk,kn);
%MRATIOPLOTS Plot forward vs backward mass flux ratios.

poben = poben*1e-5;
figure('Name','Gaseous mass flux ratio');
plot(poben,gas,'k+');
title(['Gaseous mass flux ratio, ' pstr]);
xlabel('p_1 [bar]');
ylabel('m_{forward}/m_{backward}');

figure('Name','Viscous and Knudsen mass flux ratios');
plot(poben,msk,'rd',poben,kn,'r+',poben,msv,'ks',poben,vis,'kx');
title(['Mass flux ratios, ' pstr]);
xlabel('p_1 [bar]');
ylabel('m_{forward}/m_{backward}');
legend('knudsen, mstack.','knudsen, calc.','viscous, mstack.','viscous, calc.');
legend('boxoff');
ylim([0.999 1.001]);


function computeandplot(deltap,poben,T1,s,pms,pmr)
%COMPUTEANDPLOT Compute and plot a data set.

pstr = sprintf('\\Delta p = %.1f bar',deltap*1e-5);

[fmsv,fmsk,fvis,fkn,fgas] = compute(deltap,poben,T1,s,pms);
[rmsv,rmsk,rvis,rkn,rgas] = compute(deltap,poben,T1,s,pmr);
mfluxplot(['Forward direction, ' pstr],poben,fmsv,fmsk,fvis,fkn,fgas);
mfluxplot(['Backward direction, ' pstr],poben,rmsv,rmsk,rvis,rkn,rgas);
mratioplots(pstr,poben,rgas./fgas,rmsv./fmsv,rvis./fvis,rmsk./fmsk,rkn./fkn);
