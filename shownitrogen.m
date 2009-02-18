% Data:
%[substancename datamemname exp_id T1 p1 p2 Vflow Troom] = deal(data{:});
% [T1, Troom] = C; [p1, p2] = bar; [Vflow] = ml/min; 
% Now cycle through the membranes and display nitrogen data:
sname = 'nitrogen';
lendata = size(p1,1);
s = substance(sname); f = fmodel('homogeneous');
%topology = 'porousround';
%topology = 'tube';
plots = 0;
for m = 1:size(memname,1)
  dothis = ...
      input(sprintf('Display membrane %s?  1(yes)/0 [1]: ',memname{m}));
  if isempty(dothis), dothis = true; end
  if ~dothis, continue; end

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
  lenfound = size(indices,2);
  pred = p1(indices);
  pmean = (pred + p2(indices)) / 2;
  psat = T1(indices); % allocate psat plus reuse as  T1 [K] below:
  for k = 1:lenfound
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

  Tr = Troom(sorted);
  % Set room temperature, if a value is missing
  indices = isnan(Troom(sorted)); % Reuse indices!
  if any(indices) % Some Troom-values are NaN
    if ~any(~indices) % All values are NaN!
      %disp('Room temperature set to 25 C.');
      disp(sprintf('Membrane %s: min(T1) = %0.2f C, max(T1) = %0.2f C.',...
        memname{m},min(T1(sorted))-273.15,max(T1(sorted))-273.15));
      %Tr = input(sprintf(['Set room temperature [C] for membrane %s,'...
      %	' exp. %s: '],memname{m},exp_id{sorted(1)}));
      Tr = mean(T1(sorted));
      disp(sprintf('Room temperature set to %0.2f C.',Tr-273.15));
      Tr = ones(lenfound,1) * (Tr);%+273.15);
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
  mem = membrane(poredia(m),eps(m),1.38,model{m},tau(m),beta(m),L(m));
  %kn_nu = 3*sqrt(pi/(8*s.R)) / (sqrt(T)*mem.dia);
  kn_nu = 3*sqrt(pi/(8*s.R)) / mem.dia;
  %mem2 = membrane(poredia(m),eps(m),1.38,model{m},1,beta(m),L(m));
  if beta(m) >= 8.1
    beta2 = 0.9*beta(m);
  else
    beta2 = 1.1*beta(m);
  end
  mem3 = membrane(poredia(m),eps(m),1.38,model{m},tau(m),beta2,L(m));
  flow_flux = area(m)*60*1e6; % mass flow [mg/min] over mass flux [kg/m2s]
  %L_kappa = L(m)/mem.kappa;
  % delete old calc1, preallocate array of correct length
  data = pred; calc1 = pred; idgas = pred;
  %calc2 = pred; calc3 = pred;% shouldbe1 = pred;
  for k = 1:lenfound
    ind = sorted(k);
    delp = p1(ind) - p2(ind);  vol = s.v(Tr(k),101300);
    %one_mflux = s.nug(T1(ind),pmean(k)) * L_kappa / delp;
    %disp(sprintf('p1 = %0.3f, p2 = %0.3f bar: v = %0.2f m3/kg.',...
    %  p1(ind)*1e-5,p2(ind)*1e-5,vol));
    %shouldbe1(k) = Vflow(ind) * one_mflux / (vol*flow_flux);
    %permeance(k) = Vflow(ind) * L(m) / (s.v(Tr(k),101300)*((p1(ind)-p2(ind)));
    data(k) = Vflow(ind) * L(m) / (flow_flux*vol*delp);
    %mflux = mnum(T1(ind),p1(ind),p2(ind),0,s,mem,f);
    %calc1(k) = mflux*one_mflux;
    calc1(k) = mnum(T1(ind),p1(ind),p2(ind),0,s,mem,f)*L(m)/delp;
    idgas(k) = mem.kappa * ( 1/s.nug(T1(ind),pmean(k)) ...
      + mem.beta*kn_nu/sqrt(T1(ind)) );
    %calc2(k) = mnum(T1(ind),p1(ind),p2(ind),0,s,mem2,f)*L(m)/delp;
    %calc3(k) = mnum(T1(ind),p1(ind),p2(ind),0,s,mem3,f)*L(m)/delp;
%    disp(sprintf(' %2u: %0.3f  %0.3f  %0.3f  %2.2f  %0.4g  %0.4g',...
%      k,p1(ind)*1e-5,p2(ind)*1e-5,pmean(k)*1e-5,T1(ind)-273.15,...
%      s.mug(T1(ind)),s.nug(T1(ind),pmean(k))));
  end

  % Plot
  plots = mod(plots,8);
  if plots == 0
    cfig = figure;
    set(cfig,'PaperPosition',[0.6345 0.6345 19.7150 28.408],...
      'PaperPositionMode','manual');
    %k = figure;
    %set(k,'PaperPosition',[0.5 1 get(k,'PaperSize')-[1 1.5]])%,...
    %  'DefaultAxesPosition',[0 0 1 1]);
  end
  plots = plots + 1;
  set(0,'CurrentFigure',cfig); subplot(4,2,plots);
  %plot(pred,data,'k*',pred,calc1,'k+-',pred,calc2,'k--',pred,calc3,'k:');
  plot(pred,data./idgas,'k*',pred,calc1./idgas,'k+-');
  title(sprintf('{\\bf%s}. \\tau = %0.2f; \\beta = %0.2f.',...
    memname{m},tau(m),beta(m)));
  xlabel(xlabelstr);
  %ylabel('permeance\times\mu/\kappa [-]');
  %ylabel('permeance [s]');
  ylabel('m / m_{id.gas}');
  %ylabel('Vflux [ml/min]');
  %legend('data',sprintf('calc. \\times %u',calcfac));

end
