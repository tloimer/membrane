% MEMBRANESALT Berechne tau und beta mittels linearer Regression an N2-Daten.
% Vergleiche Fehlerabschätzung (rstd) bei skalierung mit und ohne Kn/nu.
% Um ein wenig Sicherheit bei der Konstruktion der Fehlerbalken für die
% Flüssigkeitsströmung zu haben.

%filenameout = [mfilename '.tsv'];
%disp(['matlab3/' mfilename '. Schreibt: ' filenameout]);
read0905;
load results1102;
% Den Header aus membranesguessed brauchen wir auch
fid = fopen('../data/membranesguessed');
%headerout = fscanf(fid,'%s\n',1,'Delimiter','\n','Whitespace','');
headerout = textscan(fid,...
  '%s%*s%*s%s%s%s%s%s%s%s',1,'Delimiter','\t','Whitespace','');
for m = 1:length(headerout), headerout{m} = headerout{m}{1}; end
headerout = sprintf('%s\t',headerout{:});
headerout = sprintf('%srstd\trlen\trxmean\trsumsq',headerout);
fclose(fid);
% CALCTAUBETA Berechne tau und beta durch lineare Regression an N2-Daten.
% Mit Darcy'schem Gesetz, massflux = (kappa/nu_app) (del p/L), del p = p1 - p2;
% und nu_app = nu/(1 + beta Kn), Kn = lambda/d = 3nu\sqrt(pi/(8RT))/d,
% siehe Gl. (3), Loimer (2007), sowie perm = massflux*L/(delp),
%   perm = mem.kappa/(tau nu) + (beta*mem.kappa/tau)*kn_nu.
% Zuerst wurde für die Bestimmung von tau und beta lscov verwendet
%   1 = (mem.kappa/tau)/(perm*nu) + (beta*mem.kappa/tau)*(kn_nu/perm),
% mit c1c2(1) = mem.kappa/tau, c1c2(2) = beta*c1c2(1).
% Einfacher ist aber
%  1/nu = tau*(perm/mem.kappa) - beta*kn_nu.
%
% Lineare Regression verlangt besser
%  perm/(kn_nu*mem.kappa) = (1/tau)/(kn_nu*nu) + (beta/tau)
% Mit der Konstante mem.kappa, fast einer Konstante kn_nu und den bei jeder
% Messung veränderlichen perm und nu, mit nu proportional zu 1/pmean. Damit
% entspricht in der üblichen Statistik-Notation
%  hat alpha (oder hat beta_0) = (beta/tau),
%  hat beta (oder hat beta_1) = (1/tau).

%[substancename datamemname exp_id T1 p1 p2 Vflow Troom] = deal(data{:});
% [T1, Troom] = C; [p1, p2] = bar; [Vflow] = ml/min; 
%read0905;
% Liest folgende Variablen aus ../data/data0905.tsv ein:
%   L, T1, T12tc, T1tc, T2tc, Troom, Vflow, area, beta, datamemname, eps,
%   exp_id, memdia, memname, model, p1, p2, poredia, substancename, tau
%load results1102;
% enthält: kap_kc, kappac, psat1, pK1, T2calc, mexp, mgas, miso, mcalc, fl0,
%   kappa, kapKK, kapll, n, Ccc, mguess, mlinp1sat

% Now cycle through the membranes and display nitrogen data:
sname = 'nitrogen';
lendata = size(p1,1);
s = substance(sname); f = fmodel('homogeneous');
%topology = 'porousround';
topology = 'tube';
%tau = 1; beta = 8.1;
dothis = true;
lenm = size(memname,1);
daten = zeros(lenm,6);
rdaten = zeros(lenm,8);
rout = NaN(lenm,6);
for m = 1:lenm

  % Collect the nitrogen data for this membrane
  indices = []; % Preallocation is necessary.
  for k = 1:lendata
    if strcmp(memname{m},datamemname{k}) && strcmp(sname,substancename{k})
      indices = [indices k];
    end
  end

  if isempty(indices)
    disp(sprintf('No nitrogen data found for membrane %s.',memname{m}));
    continue;
  end

  % Assign pred [-] and pmean [Pa] , we probably need it anyway.
  rlen = size(indices,2);
  pred = p1(indices);
  pmean = (pred + p2(indices)) / 2;
  psat = T1(indices); % allocate psat plus reuse as  T1 [K] below:
  for k = 1:rlen
    psat(k) = s.ps(psat(k)); % ps is not vectorizable
  end
  if ~any(~isfinite(psat)) % if no infinite value is encountered in psat
    pred = pred/psat;
    xlabelstr = 'p_{red} [-]';
  else
    pred = pred*1e-5;
    xlabelstr = 'p_1 [bar]';
  end
  % else, pred stays p1 [bar].

  % Sort after pred
  % or comment out to sort after pmean.
  pred = pmean*1e-5; xlabelstr = 'p_{mean} [bar]';
  [pred sorted] = sort(pred);
  pmean = pmean(sorted);
  sorted = indices(sorted);  % Reused sorted, pmean.
  disp(sprintf('Membrane %s: min(T1) = %0.2f C, max(T1) = %0.2f C.',...
    memname{m},min(T1(sorted))-273.15,max(T1(sorted))-273.15));
  T1mean = mean(T1(sorted));

  Tr = Troom(sorted);
  % Set room temperature, if a value is missing
  indices = isnan(Troom(sorted)); % Reuse indices!
  if any(indices) % Some Troom-values are NaN
    if ~any(~indices) % All values are NaN!
      %disp('Room temperature set to 25 C.');
      %disp(sprintf('Membrane %s: min(T1) = %0.2f C, max(T1) = %0.2f C.',...
      %  memname{m},min(T1(sorted))-273.15,max(T1(sorted))-273.15));
      %Tr = input(sprintf(['Set room temperature [C] for membrane %s,'...
      %	' exp. %s: '],memname{m},exp_id{sorted(1)}));
      Tr = T1mean;
      disp(sprintf('No room temperature; Set to mean(T1), %0.1f C.',Tr-273.15));
      Tr = ones(rlen,1) * (Tr);%+273.15);
    else
      % Calculate mean
      Tr(find(indices)) = mean(Tr(find(~indices)))
    end
  end

%  disp(sprintf('%s',memname{m}));
%  disp(' #   p1     p2     pmean  T1     mu(T1)     nu(T1,pmean)');

  % And calculate for each data point
  % Because del p/L = massflux*nu/kappa, the normalized permeability should
  % ideally be one:
  %   massflux * L*nu/(del p*kappa) = 1.
  % Commonly, one sets nu = mu*R*T/pmean and plots the permeance, expecting it
  % to be linear in pmean:
  % massflux*L/(del p) = pmean * kappa/(mu*R*T)
  mem = membrane(poredia(m),eps(m),1.38,topology,1,0,L(m));
  flow_flux = area(m)*60*1e6; % mass flow [mg/min] over mass flux [kg/m2s]
  L_kappa = L(m)/mem.kappa;
  % delete old calc1, preallocate array of correct length
  perm = pred; nu = pred; kn_nu = pred;
  %data = pred; calc1 = pred; shouldbe1 = pred; calc2 = pred;
  for k = 1:rlen
    ind = sorted(k);
    delp = p1(ind) - p2(ind);  vol = s.v(Tr(k),101300);
    %one_mflux = s.nug(T1(ind),pmean(k)) * L_kappa / delp;
    %disp(sprintf('p1 = %0.3f, p2 = %0.3f bar: v = %0.2f m3/kg.',...
    %  p1(ind)*1e-5,p2(ind)*1e-5,vol));
    %shouldbe1(k) = Vflow(ind) * one_mflux / (vol*flow_flux);
    %permeance(k) = Vflow(ind) * L(m) / (s.v(Tr(k),101300)*((p1(ind)-p2(ind)));
    perm(k) = Vflow(ind) * L(m) / (flow_flux*vol*delp);
    nu(k) = s.nug(T1(ind),pmean(k));
    kn_nu(k) = 3*sqrt(pi/(8*s.R)) / (sqrt(T1(ind))*mem.dia);
    %mflux = mnum(T1(ind),p1(ind),p2(ind),0,s,mem,f);
    %calc1(k) = mflux*one_mflux;
%    calc2(k) = mnum(T1(ind),p1(ind),p2(ind),0,s,mem,f)*L(m)/delp;
%    disp(sprintf(' %2u: %0.3f  %0.3f  %0.3f  %2.2f  %0.4g  %0.4g',...
%      k,p1(ind)*1e-5,p2(ind)*1e-5,pmean(k)*1e-5,T1(ind)-273.15,...
%      s.mug(T1(ind)),s.nug(T1(ind),pmean(k))));
  end

% Calculate tau and beta by least square fit, scale mflux*L/(p1-p2) to 1.
%c1c2 = lscov([1./(nu.*perm) kn_nu./perm],ones(rlen,1));
%c1c2 = lsqnonneg([1./(nu.*perm) kn_nu./perm],ones(rlen));
%tauc = mem.kappa/c1c2(1);
%betac = c1c2(2)/c1c2(1);
disp(sprintf('Membrane %s, %u Messungen',memname{m},rlen));
%disp(sprintf('Korrelation (1/tau, beta/tau): tau = %0.4f, beta = %0.4f.',...
%  tauc,betac));

% Die einfachere Version,  1/nu = tau*(perm/mem.kappa) - beta*kn_nu.
%[tb stderr mse] = lscov([perm/mem.kappa -kn_nu],1./nu);
%disp(sprintf('Länge tb, stderr, mse: %u, %u, %u. Error mse: %.4g',...
%  length(tb),length(stderr),length(mse),mse));
%disp(sprintf(...
%  'Korrelation (tau, -beta):      tau = %.4f(%.4f), beta = %.4f(%.4f).',...
%  tb(1),stderr(1),tb(2),stderr(2)));

% Regression:
%  perm/(kn_nu*mem.kappa) = (1/tau)/(kn_nu*nu) + (beta/tau)
xi = 1./(nu.*kn_nu);
[b,bint,r,rint,stats] = ...
  regress(perm./(mem.kappa.*kn_nu),[ones(rlen,1) xi]);
meankn_nu = mean(kn_nu);
disp(['Regression,  '...
    '1/tau   beta/tau   rstd   kn_nu  (max(kn_nu)-min(kn_nu))/kn_nu']);
disp(sprintf('              %.4g  %.4g  %.4g  %.4g  %.4f%%',...
  b(2),b(1),sqrt(stats(4)),meankn_nu,(max(kn_nu)-min(kn_nu))*100/meankn_nu));
%tr = 1/b(2); trlo = 1/bint(2,2); trhi = 1/bint(2,1);
%br = b(1)*tr; brlo = bint(1,1)*trlo; brhi = bint(1,2)*trhi;
%disp(sprintf(' tau = %0.4f + %0.4f - %0.4f;    1/tau = %0.4f +/- %0.4f',...
%  tr,trhi-tr,tr-trlo, b(2),bint(2,2)-b(2)));
%disp(sprintf('beta = %0.4f + %0.4f - %0.4f; beta/tau = %0.4f +/- %0.4f',...
%  br,brhi-br,br-brlo, b(1),bint(1,2)-b(1)));
%if tr < 0
%  [b,bint,r,rint,stats] = regress(perm./(mem.kappa.*kn_nu),ones(rlen,1));
%  br = b; tr = 1;
%  disp('Erzwungen: tau = 1');
%end

xialt = 1./nu;
[ba,binta,ra,rinta,statsa] = ...
  regress(perm./mem.kappa,[meankn_nu*ones(rlen,1) xialt]);
disp(['Regr. alt.,  '...
    '1/tau   beta/tau   rstd   rstd(1)*mean(kn_nu)/rstd(2)-1']);
disp(sprintf('              %.4g  %.4g  %.4g  %.4g  %.4g',...
  ba(2),ba(1),sqrt(statsa(4)),sqrt(stats(4))*meankn_nu/sqrt(statsa(4))-1));
%disp(sprintf('Regressionalternativ: stats4 = %g, statsalt = %g',stats(4),statsa(4)));
%disp(sprintf('Regressionalternativ: xi = %g -- %g ',min(xi),max(xi)));
%disp(sprintf(' tau = %0.4f\nbeta = %0.4f',1/ba(2),ba(1)/ba(2)));

disp('Somit ist 1/tau +/- s, sigma-Daten:');

rsumsq = var(xi,1)*rlen;
rstd = sqrt(stats(4)/(var(xi,1)*rlen));
disp(sprintf('Regression: %.4g +/- %.4g (%.1f%%), +/- %.4g',...
  b(2),rstd,rstd/b(2)*100,...
    sqrt(stats(4)*(1+1/rlen))/mean(perm./(mem.kappa.*kn_nu))));
%  b(2),rstd,rstd/b(2)*100,sqrt(stats(4)*meankn_nu^2+stats(4)/rsumsq)));
rsumsq = var(xialt,1)*rlen;
rstd = sqrt(statsa(4)/(var(xialt,1)*rlen));
disp(sprintf('Regr. alt.: %.4g +/- %.4g (%.1f%%), +/- %.4g',...
  ba(2),rstd,rstd/ba(2)*100,...
    sqrt(statsa(4)*(1+1/rlen))/mean(perm./mem.kappa)));
%  ba(2),rstd,rstd/ba(2)*100,sqrt(statsa(4)+statsa(4)/rsumsq)));
% Rechnen wir besser mit korrekter Statistik
% Zuerst für den Output
%rstd = sqrt(stats(4));
%rsumsq = var(xi,1)*rlen;
%rxmean = mean(xi);
%rout(m,:) = [tr br rstd rlen rxmean rsumsq];

%% nun die verbesserte Statistik
%varb2 = stats(4)/rsumsq;
%taur = tr*(1+varb2*tr^2); betar = br*(1+2*varb2*tr^2);
%disp(sprintf('Regression (1/tau, beta/tau): tau = %0.4f, beta = %0.4f.\n',...
%  taur,betar));
%daten(m,:) = [tauc betac tr br taur betar];
%
%% Regressionsdaten
%varb_tau = stats(4)*(1/rlen + rxmean^2/rsumsq);
%xp = 1./(s.nug(298.15,2.5e5)*3*sqrt(pi/(8*s.R))/(sqrt(298.15)*mem.dia));
%varb_txp = stats(4)*(1/rlen + (xp - rxmean)^2/rsumsq);
%varmxp = stats(4) + varb_txp;
%rdaten(m,:) = [taur betar 1/tr b(1) sqrt([varb2 varb_tau varb_txp varmxp])];
%
end % for m = 1:size(memname,1)

%% Und die Zusammenfassung
%disp(...
% 'Membrane    tau     beta    tauc    betac   tauroh  betar   taur    betar');
%for m = 1:lenm
%  disp([sprintf('%-11s',memname{m}) ...
%    sprintf('%6.3f  ',tau(m),beta(m),daten(m,:))]);
%end
%disp(sprintf('\nAnalyse der Regression'));
%disp(...
%'Membrane    taur    betar   1/tau   bet/t   s(1/t)  s(b/t)  s(btxp) s(mxp)');
%for m = 1:lenm
%  disp([sprintf('%-11s',memname{m}) ...
%    sprintf('%6.3f  ',rdaten(m,:))]);
%end
%
%% Manuelle Änderungen
%
%disp(sprintf('\nÄnderungen der Membranen'));
%m = find(strcmp('11A',memname)); rout(m,logical([0 0 1 0 1 1])) = NaN;
%disp([sprintf('  %s: ',memname{m}) sprintf(' %g',rout(m,:))]);
%aendere = {'11A' 1 0.25; 'HP11A' 1 0.25; '33A' 0.94 7.32; 'HP55B' 1 7.5};
%for k = 1:size(aendere,1)
%  m = find(strcmp(aendere{k,1},memname));
%  rout(m,1) = aendere{k,2}; rout(m,2) = aendere{k,3};
%  disp(sprintf('%8s: tau = %-.3g, beta = %-.3g.', ...
%    memname{m},rout(m,1),rout(m,2)));
%end

% Write output
% Output-header
%fid = fopen(filenameout,'w');
%fprintf(fid,'%s\n',headerout);

% Schleife zu Ausgabedaten: Aus read0905
%[memname,poredia,L,memdia,eps,model,tau,beta] = deal(data{:});
%poredia = 1e-9*poredia; L = 1e-3*L; area = 1e-6*memdia.^2*pi/4; eps = 1e-2*eps;
%rout(m,:) = [tr br rstd rlen rxmean rsumsq];
%for m = 1:lenm, fprintf(fid,...
%  '%s\t%.0f\t%.3g\t%.4g\t%.4g\t%s\t%.2f\t%.2f\t%.3g\t%u\t%.3g\t%.3g\n',...
%  memname{m},poredia(m)*1e9,L(m)*1e3,memdia(m),eps(m)*100,model{m},rout(m,:));
%end
%fclose(fid);
