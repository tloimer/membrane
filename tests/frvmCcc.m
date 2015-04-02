% Massenstrom über Ccc, kappa/kappa_l = 0.265, 0.9 < P1 < 1.0
read1103;
% Liest folgende Variablen aus ../data/data0905.tsv ein:
%   L, T1, T12tc, T1tc, T2tc, Troom, Vflow, area, beta, datamemname, eps,
%   exp_id, memdia, memname, model, p1, p2, poredia, substancename, tau
load /home/tloimer/projects/11jms/matlab3/results1105
% enthält: kap_kc, kappac, psat1, pK1, T2calc, mexp, mgas, miso, mcalc, fl0,
%   kappa, kapKK, kapll, n, Ccc, mguess, mlinp1sat

global VERBOSE
figure('Name',mfilename);

% FigureDefaults
set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',8,...
   'DefaultTextFontSize',8);
set(0,'DefaultFigurePaperUnits','points',...
  'DefaultFigurePaperPositionMode','manual',...
  'DefaultFigurePaperPosition',[0 0 262 189]); % sicherer: 240 points
% 248 186 - gibt 230-235 x 180-183
% soll <= 251 sein
set(0,'DefaultLineColor','k','DefaultLineMarkerEdgeColor','k');
  %,'DefaultLineMarkerFaceColor','k');

% Hier extrahieren wir zuerst die Daten, dann suchen wir die
% Membraneigenschaften heraus und berechnen damit die theoretischen Linien.
% In fmvskkl.m geschah dies umgekehrt. Fürs Daten filtern müssen wir die
% Temperaturen wissen.

% Temperaturen

% Einkopiert aus tempdifferenz.m
isDT12meas = isfinite(T12tc);
isT1T2meas = isfinite(T1tc) & isfinite(T2tc);

% allocate array T12meas; assign T12meas from two sources
lendata = size(p1,1);
T12meas = NaN(lendata,1);
for i = 1:lendata
  if isDT12meas(i)
    T12meas(i) = T12tc(i)/(-41);
  elseif isT1T2meas(i)
    T12meas(i) = (T2tc(i)-T1tc(i))/41;
  end
end

T12calc = T1 - T2calc;
isT12ok = T12meas < 1.31*T12calc & T12meas > 0.3*T12calc;
clear('isDT12meas','isT1T2meas');
% Ende Kopie aus tempdifferenz.m.

% Daten aussuchen
isvapor = ~strcmp(substancename,'nitrogen');
P1 = (p1-pK1+n.*(p1-p2))./(psat1-pK1+n.*(p1-p2));
kappall = kappa./kapll; % reuse kappall
isexp = isvapor & isT12ok & kappall>0.265 & kappall<0.266 & P1>0.9;

% Membraneigenschaften
thismem = unique(datamemname(isexp));
thissubstance = unique(substancename(isexp));

if length(thissubstance) > 1, error; end

i = find(strcmp(thismem,memname));
if ~isscalar(i), error('More than one membrane'); end

mem = membrane(poredia(i),eps(i),1.38,model{i},tau(i),beta(i),L(i));

f = fmodel('homogeneous');
s = substance(thissubstance{1});
theta=0;

% Abkürzungen
meg = mexp(isexp)./mgas(isexp);

lendata = sum(isexp);
ind = find(isexp);
mstacktest = zeros(lendata,1); % All other variables are column vectors!
for i = 1:lendata
  ii = ind(i);
  if VERBOSE > 0
    fprintf('%s: # %u, altes mnum: mguess %.3g, mcalc %.3g.\n',...
      upper(mfilename),ii,mguess(ii),mcalc(ii));
  end
  mstacktest(i) = mstack(T1(ii),p1(ii),p2(ii),theta,s,mem,f);
%  mstacktest(i) = mnum(T1(ii),p1(ii),p2(ii),theta,s,mem,f);
end

hl = plot(Ccc(isexp),mcalc(isexp),'ok',Ccc(isexp),mstacktest,'+k');

xlim([0 10]);
%ylim([0 10]);
xlabel('C_{cc}');
ylabel('m','Rotation',0);
legend('calc. mnum','calc. mstack','Location','NorthWest');

legend('boxoff');
box('on');

print('-deps2',[mfilename '.eps']);
