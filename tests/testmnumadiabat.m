%TESTMNUMADIABAT Test asym.m, mstackstruct.m and mnumadiabat

global VERBOSE;
VERBOSE = 0;
% Verbositiy levels:
% 0 - nothing
% 1 - print diagnostic output from internally available data
% 2 - iteration information (fzero(?), findzero)
% 3 - re-calculate or calculate additional data;
%     calculate and plot residual pressure vs. mass flux in findinterval.m

% To force this file to be called from the directory tests/, because the files
% flow12vsasym.m and flow12andasym.m it calls also are within tests/.
%[~,name] = fileparts(cd);
%if name ~= 'tests', error('Run from directory tests/.'); end
%clear('name');

if ~exist('substance.m'), addpath('../program'); end

% A few membranes
mem1 = membrane(10e-9,0.5,1.2,'porousround',8.1,1,1e-3);
mem2 = membrane(50e-9,0.5,0.8,'porousround',8.1,1,0.5e-3);
mem3 = membrane(300e-9,0.5,0.8,'porousround',8.1,1,2e-3);
%mem3 = membrane(1000e-9,0.95,10,'porousround',8.1,1,2e-3);
mem4 = membrane(600e-9,0.5,0.8,'porousround',8.1,1,2e-3);

but = substance('butane');
iso = substance('isobutane');
n2 = substance('nitrogen');

f = fmodel('plug');

fprintf(['Large pores, two-phase flow. This triggered an error. Findzero was used\n'...
'(commit e12964f), solver accuracy switching was removed (commit b1b062), nothing\n'...
'helped. Finally, a signal was raised if a pure flow integration terminated\n',...
'prematurely (commit 6652806). Then, e.g., p should be equal to pk(T), but was\n'...
'only to within computer accuracy. That signal helped.\n']);

fprintf(['Here come the two cases where iteration initially failed, due to wrong\n'...
'integration: 2ph-flow was not commenced after vapor flow terminated prematurely.\n'...
'The mass flux is set, therefore a large residual pressure remains.\n']);

% Set up the membranestruct
mstmp = mstackstruct(0,mem3,f);

% Copy from mnumadiabat
% Calculate the downstream temperature
T1 = 290.5; p1 = 1.9e5; p2 = 1.0e5;
T2 = but.intjt(T1,p1,p2);
% Set the downstream state, downstreamstate(T2,p2,a2,q2,s,m)
state2 = downstreamstate(T2,p2,[],0,but); % thus, m is empty

% Copy some values to the membrane struct;
mstmp.T1 = T1;
mstmp.p1in = p1;
mstmp.T2 = T2;
mstmp.p2 = p2;
mstmp.a2 = state2.a;
mstmp.q2 = state2.q; % = 0

% Calculate and set flowsetup-structures
mstmp = mstmp.writeflowsetups(T1,T2,but,mstmp);
solvacc = solverstruct('accurate');
solvacc.writesolution = true;
solvacc.fullsolution = true;

msok = mstmp;
msok.m = 0.0017;
[p1ok,msok] = asym(msok.m,state2,msok,solvacc);
msok.p1sol = p1ok;
mstmp.printsolution(msok);

msirr = mstmp;
msirr.m = 0.00246995;
[p1irr,msirr] = asym(msirr.m,state2,msirr,solvacc);
msok.p1sol = p1irr;
mstmp.printsolution(msirr);

fprintf('To resume, type "return".\n');
%keyboard;

fprintf(['\nThese iterations failed, crude and accurate solver,'...
'for pore diameter\n300 nm and 600 nm, respectively.\n'...
'pres vs mflux from findinterval, further below, showed spikes.\n']);

VERBOSE = 3;
for dia = [300 600]
  memtmp = membrane(dia*1e-9,0.5,0.8,'porousround',8.1,1,2e-3);
  fprintf('\nCrude solver in asym, pore dia = %d nm\n', dia);
  flow12vsasym(290.5,1.9e5,1.0e5,0,but,memtmp,f,'crude')
  fprintf('\nAccurate solver in asym\n');
  flow12vsasym(290.5,1.9e5,1.0e5,0,but,memtmp,f)
  VERBOSE = 0;
end

fprintf('To remove figures, type "close all". To resume, type "return".\n');
keyboard;
clear memtmp mstmp msok msirr T1 T2 p1 p2 state2 solvacc p1ok p1irr dia;

fprintf(['\nNow crude vs. accurate solver shows little difference,\n',...
'not large differences as before.\n']);

dia = [200:20:600];
% allocate
mflow = dia; mcrude = dia; maccurate = dia;
pflow = dia; pcrude = dia; pacc = dia;
for i = 1:length(dia)
  memtmp = membrane(dia(i)*1e-9,0.5,0.8,'porousround',8.1,1,2e-3);
  [mflow(i),mcrude(i),maccurate(i),pflow(i),pcrude(i),pacc(i)] ...
      = flow12andasym(290.5,1.9e5,1e5,0,but,memtmp,f);
  fprintf('%d nm: flow12 = %.4g g/m2s, crude %.1f%%, accurate %.1f%%\n',...
	  dia(i),mflow(i)*1e3,100*mcrude(i)/mflow(i),100*maccurate(i)/mflow(i));
  fprintf('p1calc - p1: %+.0f Pa,   %+.0f Pa,   %+.0f Pa\n',...
	  pflow(i),pcrude(i),pacc(i));
end
figure('Name','Massflux - diameter');
plot(dia,mflow*1e3,'ks',dia,mcrude*1e3,'ro',dia,maccurate*1e3,'k+');
xlabel('pore diameter [nm]');
ylabel('mass flux [g/m2s]');
xlim([190 610]);
legend('flow12','crude','accurate');
legend('Location','NorthWest');
legend('boxoff');
figure('Name','Residual pressure');
plot(dia,pflow,'ks',dia,pcrude,'ro',dia,pacc,'k+');
xlabel('pore diameter [nm]');
ylabel('p_{1,calc} - p_1 [Pa]');
xlim([190 610]);
legend('flow12','crude','accurate');
legend('Location','Best');
legend('boxoff');

fprintf('To remove figures, type "close all". To resume, type "return".\n');
keyboard

fprintf('\nThe (nearly) systematic testing of all possible cases starts.\n\n');

fprintf('\nA liquid film close to saturation\n');
flow12vsasym(290.5,1.9e5,1e5,80,but,mem1,f);

fprintf('\nA liquid film at saturation\n');
% Had to allow minute negative q1 in asym>heatfluxcriterion
flow12vsasym(290.5,but.ps(290.5),1e5,80,but,mem1,f);

fprintf('\nStill a film, ideally wetting\n');
% Found an error in asym>interface_vapliq, temperature was 42 K too high!
% corrected fs.hvapKraw(T2,p2,...) fs.hvapKraw(T2,p1,...);
% also tested to revert last change (heatfluxcriterion) - both are necessary
flow12vsasym(290.5,1.9e5,1e5,0,but,mem1,f)

fprintf('To remove figures, type "close all". To resume, type "return".\n');
keyboard

fprintf('\nLarger pore diameter, no liquid film any more\n');
flow12vsasym(290.5,1.9e5,1e5,0,but,mem2,f)

fprintf('\nAgain small pores, a liquid film and liquid in the entire membrane.\n');
flow12vsasym(290.5,1.9e5,1.7e5,0,but,mem1,f)

fprintf('\nLarger pores, no film, liquid in the entire membrane.\n');
flow12vsasym(290.5,1.9e5,1.82e5,0,but,mem2,f)

fprintf('To remove figures, type "close all". To resume, type "return".\n');
keyboard

fprintf('\nLarge pores, two-phase flow\n');
flow12vsasym(290.5,1.9e5,1.0e5,0,but,mem3,f);

fprintf(['\nIncreasing the downstream pressure,the two-phase flow region becomes\n'...
'longer and longer. At p2 = 188.6 kPa two-phase flow is in nearly the entire\n'...
'membrane, at p2 = 188.7 there is liquid in the entire membrane.\n']);

fprintf('\nTwo-phase flow in nearly the entire membrane\n');
flow12vsasym(290.5,1.9e5,1.886e5,0,but,mem3,f)

fprintf('\nLiquid in the entire membrane\n');
flow12vsasym(290.5,1.9e5,1.887e5,0,but,mem3,f)

fprintf('To remove figures, type "close all". To resume, type "return".\n');
keyboard

fprintf('\nNon-wetting cases\n\n');

fprintf('\nAnother case that triggered an error during iteration: theta = 90\n');
fprintf('Liquid film, liquid flow through part of the membrane\n');
fprintf(['Needs, during iteration, the pressure set in state.p for '...
'two-phase flow.\nDone in commit a71d367\n']);
flow12vsasym(290.5,1.9e5,1.0e5,90,but,mem3,f);

fprintf('\nLarge pores, nevertheless a liquid film in front of the membrane\n');
flow12vsasym(290.5,1.9e5,1.0e5,120,but,mem3,f);

fprintf('\nContact angle 95 degrees, liquid flow through part of the membrane\n');
flow12vsasym(290.5,1.9e5,0.5e5,95,but,mem3,f);

fprintf('\nStill 95 degree contact angle, a liquid film and two-phase flow\n');
flow12vsasym(290.5,1.9e5,1.0e5,95,but,mem4,f);

fprintf('To remove figures, type "close all". To resume, type "return".\n');
keyboard

fprintf('\nContact angle theta = 90, saturated vapor\n');

psat1 = but.ps(290.5);

fprintf('\nLiquid film, liquid flow through part of the membrane\n');
flow12vsasym(290.5,psat1,1e5,90,but,mem2,f);

fprintf('\nStill a haze of a liquid film\n');
flow12vsasym(290.5,psat1,1e5,90,but,mem3,f);

fprintf('\nTwo-phase flow, runs after commit 80513d8\n');
flow12vsasym(290.5,psat1,1.0e5,90,but,mem4,f);

fprintf('To remove figures, type "close all". To resume, type "return".\n');
keyboard

fprintf(['\nThe same as above, now with one ore more support layers and ',...
'additional membranes.\nRunning after commit 6047925, but temperatures far off.\n']);

mem5 = membrane(5e-6,0.8,0.8,'porousround',8.1,1,2e-3);
mem6 = membrane(50e-6,0.8,0.8,'porousround',8.1,1,2e-3);

ms1 = {{mem1 mem3} mem5};
ms2 = {{mem2 mem3} {mem5}};
ms3 = {{mem3 mem5} {mem5 mem6}};
ms4 = {{mem4 mem6} {mem6}};

fprintf('\nWetting cases\n');

fprintf(['\nFirst, the low-pressure difference cases, ',...
'with liquid all through the entire membrane\n']);

fprintf('\nLiquid film - liquid flow\n');
ms = mstackstruct(0,ms1,f);
[~,ms] = mnumadiabat(290.5,psat1,1.7e5,but,ms); ms.printsolution(ms);

fprintf('\nLiquid flow through the entire membrane\n');
ms = mstackstruct(0,ms2,f);
[~,ms] = mnumadiabat(290.5,psat1,1.82e5,but,ms); ms.printsolution(ms);

fprintf('\nLarge pressure difference cases, evaporation within the membrane\n');

fprintf('\nLiquid film - liquid flow - vapor flow\n');
ms = mstackstruct(0,ms1,f);
[~,ms] = mnumadiabat(290.5,psat1,1e5,but,ms); ms.printsolution(ms);

fprintf('\nLiquid flow - vapor flow\n');
ms = mstackstruct(0,ms2,f);
[~,ms] = mnumadiabat(290.5,psat1,1e5,but,ms); ms.printsolution(ms);

fprintf('\nTwo-phase flow - vapor flow\n');
ms = mstackstruct(0,ms3,f);
[~,ms] = mnumadiabat(290.5,psat1,1e5,but,ms); ms.printsolution(ms);
