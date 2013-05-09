function pTplot(name,fl,ispgfplot)
%PTPLOT     Plot p-T diagram.
%  PTPLOT(NAME,FLOWSTRUCT,ISPGFPLOT) plots a p-T diagram. If ISPGFPLOT is true,
%  a pgfplot-file with the name NAME.pgfplot is created.

global VERBOSE;

%fl.info.T1, .p1, .p2, .substance, .membrane, .fmodel, .flsetup
%fl.sol.T2, .T1, .len, .colors
%fl.calc.psat1, .pK1, .n

% The number of line-parts.
nflow = -fl.sol.len;

% Flattened distributions.
T = [fl.flow(1:end).T];
p = [fl.flow(1:end).p];

% The temperature range
margin = 2; step = 2;
Tmax = ceil( (max(T)+margin)/step ) * step;
Tmin = floor( (min(T)-margin)/step ) * step; 

% The pressure range [Pa]
margin = 1e4; step = 2e4;
% This is pmax if p1 <= psat(T1)
%pmax = ceil( (fl.info.substance.ps(Tmax)+margin)/step ) * pstep;
% But be robust, e.g., for non-wetting systems.
pmax = ceil( (max([p fl.info.substance.ps(Tmax)])+margin)/step ) * step;
pmin = floor( (min(p)-margin)/step ) * step;
pmax = pmax / 1e5;
pmin = pmin / 1e5;

% The saturation-pressure, pK and pK-pcap lines.
nTsat = 5;
step = (Tmax-Tmin)/(nTsat-1); % re-defined step
% Add 0.1*Tstep, to prevent failing with rounding errors.
Tsat = [Tmin:step:Tmax+0.1*step];
% Memory allocation
psat = Tsat;
pK = Tsat;
pKcap = Tsat;
% psat, pK and pKcap are written in units of bar!
for i = 1:nTsat
  psat(i) = fl.info.substance.ps(Tsat(i))/1e5;
  pK(i) = psat(i)*fl.info.flsetup.pkps(Tsat(i));
  pKcap(i) = pK(i) - fl.info.flsetup.curv*fl.info.substance.sigma(Tsat(i))/1e5;
end

if VERBOSE > 0
  fprintf('pKcap: [%.3g %.3g] bar\n', min(pKcap), max(pKcap));
end

%%%	THE PLOT	%%%
name = [name '.pgfplot'];
% pgfplots
if ispgfplot
  rangestr = sprintf(' xmin =%d, xmax = %d, ymin = %d, ymax = %d\n',...
	Tmin, Tmax, pmin, pmax);
  pid = beginpgfplot(name, ['xlabel = {$T$ [K]}, ylabel = {$p$ [bar]},\n' rangestr ...
    ' legend style={at={(0.97,0.07)},anchor=south east,cells={anchor=west}},\n'...
    ' y label style = {rotate=-90,xshift=-10bp}, width=8cm, height=6cm']);

  addcoords(pid,Tsat',psat','black, solid, thick');
  addcoords(pid,Tsat',pK','black, dashed, thin');
  addcoords(pid,Tsat',pKcap','black, dot-dashed, thin');

  for i = 1:nflow
    % to be mended
    addcoords(pid,fl.flow(i).T',fl.flow(i).p'/1e5,...
	'mark=*,mark options={scale=0.6},solid,thick');
  end
  endpgfplot(pid,...
	'$p_\mathrm{sat}$, $p_\mathrm K$, $p_\mathrm K - p_\mathrm{cap}');
else
% matlab-plots
  figure('Name',['p-T diagram:  ' name]);
  plot(Tsat,psat,'LineWidth',2,'LineStyle','-','Color',[0.8 0.8 0.8]);
  xlim([Tmin Tmax]);
  ylim([pmin pmax]);
  xlabel('T [K]');
  ylabel('p [bar]');
  line(Tsat,pK,'LineStyle','--','Color','k');
  line(Tsat,pKcap,'LineStyle','-.','Color','k');
  % hold on; plot..; hold off ginge auch, ohne xlim, ylim setzen zu müssen
  for i = 1:nflow
    line(fl.flow(i).T,fl.flow(i).p/1e5,'Color',fl.flow(i).color,...
	'Marker','*','LineStyle','-');
  end
end
