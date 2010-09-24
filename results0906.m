% readdata füllt folgende Arrays:
% [memname,poredia,L,memdia,eps,model,tau,beta] = deal(data{:});
% [substancename datamemname exp_id T1 p1 p2 Vflow Troom T1tc T2tc T12tc]...
% Temperaturleitfähigkeit hard-coded = 1.38 W/m2K.

% Arrays allokieren:
kap_kc = zeros(size(p1));
[kappa kappac kapKK kapll psat1 pK1 T2calc n Ccc mlinp1sat ...
  mexp mgas miso mcalc mguess] = deal(kap_kc);
% kap_kc = kappa/kappa_c;
% kappaK = kappaK(kappac), nicht genau kappaK! Schreit bei mehr als 1% Differenz
% mexp = experimentell gemessener Massenfluss,
% mgas = Massenfluss, falls der Dampf ein idealas Gas wäre,
% miso = isotherm, m90 = Kontaktwinkel von 90° (WS),
% mcalc = 0° Kontaktwinkel,
% mjt0 = (T1 = T2) aber nicht isotherm (Joule-Thomson Koeffizient = 0);
% unphysikalisch, wegen q2 = 0, nicht gerechnet.

% Um die Raumtemperatur zuzuschreiben, die manchmal fehlt, werden die Daten nach
% der Reihe durchgegangen, Daten mit gleicher Membrane und gleicher Substanz
% zusammengefasst.

fehler0 = [];
fehler90 = [];
f = fmodel('homogeneous');
Tnorm = 25 + 273.15;
first = 1;
lendata = size(p1,1);
while first <= lendata
  mname = datamemname{first};
  sname = substancename{first};
  next = first + 1; % first muss nicht mit sich selbst verglichen werden
  while next <= lendata && ...
    strcmp(mname,datamemname{next}) && strcmp(sname,substancename{next})
    next = next + 1;
  end
  last = next - 1;
  % Nun haben wir zwischen [first:last] gleiche Substanz und gleiche Membran.
  % Raumtemperatur berechnen, falls nötig
  noTr = isnan(Troom(first:last));
  if any(noTr)
    % mind. 1 Datenpunkt hat keine Raumtemperatur
    if all(noTr)
      % kein einziger Punkt hat eine Raumtemperatur; Wir setzen 25°C;
      Troom(first:last) = Tnorm;
    else
      % einzelne Punkte ohne Raumtemperatur wird der Durchschnitt des restlichen
      % Experimentes zugewiesen.
      Troom(first-1+find(noTr)) = mean(Troom(first-1+find(~noTr)));
    end
  end
  % Ende Raumtemperatur berechnen / zuweisen.

  % Der Spass beginnt. 
  s = substance(sname);
  ismem = strcmp(mname,memname); % ein logisches Array zur Indizierung;
  % doppelt vorkommende Membranen sollten einen Fehler geben.
  if sum(ismem) > 1, error(['Mehrere gleiche Membranen ' mname '!']); end
  mem = membrane(poredia(ismem),eps(ismem),1.38,model{ismem},...
    tau(ismem),beta(ismem),L(ismem));
  flow_flux = area(ismem)*60*1e6; % mass flow [mg/min] over mass flux [kg/m2s]
  betakn_nu = mem.beta * 3*sqrt(pi/(8*s.R))/mem.dia;
  kap_L = mem.kappa/mem.L;
  for i = first:last
    % gemessener Massenstrom
    mexp(i) = Vflow(i)/(flow_flux*s.v(Troom(i),101300));
    % berechnet, numerisch, theta = 0°.
    try
      [mcalc(i) fl] = mnum(T1(i),p1(i),p2(i),0,s,mem,f); 
      T2calc(i) = fl.sol.T2;
      fl0{i} = fl;
      mguess(i) = fl.calc.mguess;
      % f_handles müssen ohnehin die Funktionen finden; Vars. verkleinern.
      fl.info = rmfield(fl.info,{'substance' 'membrane' 'fmodel'});
      fl0{i} = fl;
      % fl0{i}.info.substance = sname; fl0{i}.info.membrane = mname;
      % fl0{i}.info.fmodel = fl0{i}.info.fmodel.name;
    catch ME
      disp(sprintf('Fehler bei Punkt %u: mcalc',i));
      fehler0 = [fehler0 i];
      for k = 1:length(ME.stack)
        disp(sprintf(' Zeile %u in %s.m',ME.stack(k).line,ME.stack(k).name));
      end
    end
    % Gasströmung, isotherm
    corrKn = betakn_nu/sqrt(T1(i));
    mgas(i) = kap_L*(p1(i)-p2(i))*( 1/s.nug(T1(i),(p1(i)+p2(i))/2) + corrKn );
    % Kondensation ist nur bei Dämpfen möglich
    if ~strcmp(sname,'nitrogen') && size(fl0,2)==i
      [psat1(i) tmp] = s.ps(T1(i)); %tmp ist dps/dT
      kappac(i) = s.nul(T1(i))*f.kmliq(mem.epsilon,mem.km,s.kl(T1(i)))...
        /(tmp*s.hvap(T1(i)));
      % kleine Kontrolle kappac(i) =...
      epsilon = 2^(-51);
      if kappac(i) ~= fl0{i}.calc.kappac && ...
        1-kappac(i)/fl0{i}.calc.kappac > epsilon
        disp(sprintf('Punkt %u: kappac ~= fl.calc.kappac, 1-kc/fl.kc = %g',...
	  i,1-kappac(i)/fl0{i}.calc.kappac));
      end
      % kleine Kontrolle mgas(i) =...
      if mgas(i) ~=fl0{i}.calc.mgas && 1-mgas(i)/fl0{i}.calc.mgas >epsilon
	disp(sprintf('Punkt %u: mgas ~=fl.calc.mgas, 1-mgas/fl.mg = %g',...
	  i, 1-mgas(i)/fl0{i}.calc.mgas));
      end
      kap_kc(i) = mem.kappa/kappac(i);
      kappa(i) = mem.kappa;
      kapKK(i) = fl0{i}.calc.kapKK;
      kapll(i) = fl0{i}.calc.kapll;
      Ccc(i) = fl0{i}.calc.Ccc;
      n(i) = fl0{i}.calc.n;
      mlinp1sat(i) = fl0{i}.calc.mlinp1sat;
      % isotherm, eventuell mit Kondensation
      pcap = mem.fcurv(1)*s.sigma(T1(i));
      pK1(i) = psat1(i)*exp(-pcap/(s.R*s.rho(T1(i))*T1(i)));
      % kleine Kontrolle pK1(i) =...
      if pK1(i) ~= fl0{i}.calc.pK1 && 1-pK1(i)/fl0{i}.calc.pK1 > epsilon
        disp(sprintf('Punkt %u: pK1 ~= fl.calc.pK1, 1-pK1/calc.pK1 = %g',...
	  i,1-pK1(i)/fl0{i}.calc.pK1));
        if psat1(i) ~= fl0{i}.calc.psat1
          disp(sprintf(...
	    'Punkt %u: psat1 ~= fl.calc.psat1, 1-psat1/calc.psat1 = %g',...
	    i,1-psat1(i)/fl0{i}.calc.psat1));
        end
      end
      if p1(i) >= pK1(i)
        pv_nu = (pK1(i)-p2(i))*( 1/s.nug(T1(i),(pK1(i)+p2(i))/2) + corrKn );
        pl_nu = (p1(i)-pK1(i)+pcap*(p1(i)-pK1(i))/(psat1(i)-pK1(i)))...
	  /s.nul(T1(i));
        miso(i) = kap_L*(pv_nu+pl_nu);
      else
        miso(i) = mgas(i);
      end
    else % nitrogen
      miso(i) = mgas(i);
    end %end if ~nitrogen
  end
  first = next;
end

if ~isempty(fehler0)
  disp(sprintf('Fehler für theta = 0 in %u Datenpunkten.',size(fehler0,2)));
end
if ~isempty(fehler90)
  disp(sprintf('Fehler für theta = 90° in %u Datenpunkten.',size(fehler90,2)));
end

%file = 'results0905';
dlmwrite([mfilename '.tsv'],...
  [kap_kc kappac psat1 pK1 T2calc mexp mgas miso mcalc],'\t');
save(mfilename,'kap_kc','kappac','psat1','pK1','T2calc','mexp','mgas','miso',...
  'mcalc','fl0','kappa','kapKK','kapll','n','Ccc','mguess','mlinp1sat');
