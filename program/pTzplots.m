function pTzplots(name,fl,ispgfplot)
%PTZPLOTS   Plot p-z and T-z diagrams.
%  PTZPLOTS(NAME,FLOWSTRUCT,ISPGFPLOT) plots a p-z diagram and a T-z diagram. If
%  ISPGFPLOT is true, two pgfplot-files are created with the names
%  pNAME.pgfplot and TNAME.pgfplot, respectively.

%fl.info.T1, .p1, .p2, .substance, .membrane, .fmodel, .flsetup
%fl.sol.T2, .T1, .len, .colors
%fl.calc.psat1, .pK1, .n

% The number of line-parts.
nflow = -fl.sol.len;
z3 = fl.flow(nflow).z(1);
L = fl.info.membrane.L;
zscale = fl.sol.zscale;
% Aus flow12.m
%    z3 = FL.flow(-FL.sol.len).z(1),  zscale = FL.sol.zscale,
%    z = z3 + zscale * log( (Fl.flow(-FL.sol.len).z-z3)/zscale + 1 ).
%z = [fl.flow(1:end-1).z z3+zscale*log((fl.flow(nflow).z-z3)/zscale+1)]/L;
zfront = z3 + zscale*log((fl.flow(nflow).z-z3)/zscale+1);
zflat = [fl.flow(1:end-1).z zfront] / L;

% Flattened distributions.
T = [fl.flow(1:end).T];
p = [fl.flow(1:end).p];

% The temperature range
margin = 2; step = 2;
Tmax = ceil( (max(T)+margin)/step ) * step;
Tmin = floor( (min(T)-margin)/step ) * step; 

% The pressure range [Pa]
margin = 1e4; step = 2e4;
pmax = ceil( (max(p)+margin)/step ) * step;
pmin = floor( (min(p)-margin)/step ) * step;
pmax = pmax / 1e5;
pmin = pmin / 1e5;

% z-range
step = 0.2; margin = 0.1;
Tzmin = floor((-3*zscale/L-margin)/step) * step;
pzmin = -0.1;
zmax = 1.1;

%%%	THE PLOTS	%%%
Tname = ['T' name '.pgfplot'];
pname = ['p' name '.pgfplot'];
% pgfplots
if ispgfplot
  % T-z diagram
  rangestr = sprintf(' xmin =%.1f, xmax = %.1f, ymin = %d, ymax = %d\n',...
	Tzmin, zmax, Tmin, Tmax);
  Tid = beginpgfplot(Tname,['xlabel = {$z/L$}, ylabel = {$T$ [K]},\n' rangestr ...
    ' legend style={at={(0.97,0.07)},anchor=south east,cells={anchor=west}},\n'...
    ' y label style = {rotate=-90,xshift=-10bp}, width=8cm, height=6cm']);

  addcoords(Tid,zflat',T','mark=*,mark options={scale=0.6},solid,thick');
  endpgfplot(Tid);

  % p-z diagram
  rangestr = sprintf(' xmin =%.1f, xmax = %.1, ymin = %.1f, ymax = %.1fn',...
	pzmin, zmax, pmin, pmax);
  pid = beginpgfplot(['p' name '.pgfplot'], [...
    'xlabel={$z/L$}, ylabel = {$p$ [bar]},\n' rangestr ...
    ' legend style={at={(0.97,0.07)},anchor=south east,cells={anchor=west}},\n'...
    ' y label style = {rotate=-90,xshift=-10bp}, width=8cm, height=6cm']);

  for i = 1:nflow
    addcoords(pid,fl.flow(i).z'/L,fl.flow(i).p'/1e5,...
	'mark=*,mark options={scale=0.6},solid,thick');
  end
  endpgfplot(pid);
else
% matlab-plots
  % T-z diagram
  figure('Name',['T-z diagram:  ' Tname]);
  %plot(zflat,T,'k*');
  plot(fl.flow(1).z/L,fl.flow(1).T,'Color',fl.flow(1).color,...
	'LineStyle','-','Marker','*');
  xlim([Tzmin zmax]);
  xlabel('z/L');
  ylabel('T [K]');
  for i = 2:nflow
    line(fl.flow(i).z/L,fl.flow(i).T,'Color',fl.flow(i).color,...
	'Marker','*','LineStyle','-');
  end
  % p-z diagram
  figure('Name',['p-z diagram:  ' pname]);
  plot(fl.flow(1).z/L,fl.flow(1).p/1e5,'Color',fl.flow(1).color,...
	'LineStyle','-','Marker','*');
  xlim([pzmin zmax]);
  ylim([pmin pmax]);
  xlabel('z/L');
  ylabel('p [bar]');
  % hold on; plot..; hold off ginge auch, ohne xlim, ylim setzen zu müssen
  for i = 2:nflow
    line(fl.flow(i).z/L,fl.flow(i).p/1e5,'Color',fl.flow(i).color,...
	'Marker','*','LineStyle','-');
  end
end
