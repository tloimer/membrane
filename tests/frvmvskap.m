% Massenstrom über kappa/kappa_l, Ccc > 1, 43.8 < (p1-p2)/kPa < 49.7
% 0.927 < P1 < 0.983 für k/k_l < 1;

% Bezüglich der Abschätzung von 1/tau und beta/tau kann einem gewissen
% Massenstrom nicht ein gewisser Wert von tau und beta zugewiesen werden. Es
% wird aber der Bereich von 1/tau festgestellt, der zu 95% den wahren Wert
% einschließt. Damit werden horizontale Linien als Fehlerbalken gezeichnet.
%read1103;
% Daten für jeden Datenpunkt:
% substancename datamemname exp_id T1 p1 p2 Vflow Troom T1tc T2tc T12tc
% Membrandaten;
% memname poredia L memdia eps model tau beta rstd rlen rxmean rsumsq
%load /home/tloimer/projects/11jms/matlab3/results1105
% enthält: kap_kc, kappac, psat1, pK1, T2calc, mexp, mgas, miso, mcalc, fl0,
%   kappa, kapKK, kapll, n, Ccc, mguess, mlinp1sat

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

% Die experimentellen Daten. Zuerst benötigen wir die Temperaturen, dann
% werden die weiteren Bedingungen angewendet.

% Temperaturen überprüfen

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

% extract experimental data
p12lo = 43799; p12hi = 49701;
isvapor = ~strcmp(substancename,'nitrogen');
P1 = (p1-pK1+n.*(p1-p2))./(psat1-pK1+n.*(p1-p2));
kappall = kappa./kapll; % reuse kappall
isexp = isvapor & isT12ok & ((P1>0.927 & P1<0.983 & kappall<1)|kappall>1) &...
   Ccc>1 & p1-p2>p12lo & p1-p2<p12hi;
%isexp(238)=true; % kappa/kappa_ll = 5.879, p1-p2 = 53 200 Pa.
zpunkt=238;
% Nun wird der Zusatzpunkt ununterscheidbar zu den anderen Werten
% dazugeschrieben.
isexp(zpunkt)=true;

% Calculate with mstack.
theta = 0;
f = fmodel('homogeneous');
lendata = sum(isexp); ind = find(isexp);
% Allocate vector mstacktest. Make it a column vector, like all other variables.
mstacktest = zeros(lendata,1);
for i = 1:lendata
  ii = ind(i);
  mm = find(strcmp(datamemname(ii),memname));
  if ~isscalar(mm), error([upper(mfilename) ': More than one membrane!']); end
  ms = mstackstruct(theta,  membrane(poredia(mm),eps(mm),1.38,...
				model{mm},tau(mm),beta(mm),L(mm)),  f);
  mstacktest(i) = mnumadiabat(T1(ii),p1(ii),p2(ii), ...
			      substance(substancename{ii}),ms);
end

mcg = mcalc(isexp)./mgas(isexp);

semilogx(kappa(isexp)./kapll(isexp),mcg,'ok',...
  kappa(isexp)./kapll(isexp),mstacktest./mgas(isexp),'+k');

%set(gca,'Xscale','log');
xlim([0.1 10]);
ylim([0 12]);
xlabel('\kappa/\kappa_l');
ylabel('massflux/massflux_{gas}','Rotation',90);
set(gca,'XTickLabel',[0.1 1 10]); %,...
legend('calc. 1105', 'calc. current', 'Location','NorthEast');
legend('boxoff');
%box('on');

print('-deps2',[mfilename '.eps']);
