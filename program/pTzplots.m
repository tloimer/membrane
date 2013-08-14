function pTzplots(fl,name,ispgfplot,isqplot,vpos)
%PTZPLOTS   Plot p-z, T-z and, optionally, q-z diagrams.
%  PTZPLOTS(FLOWSTRUCT) Plot a p-z and a T-z diagramm.
%  PTZPLOTS(FLOWSTRUCT,NAME,ISPGFPLOT,ISQPLOT,VPOS) plots a p-z, a T-z and,
%  optionally, a q-z diagram. NAME, ISPGFPLOT, ISQPLOT and VPOS are all optional.
%  If ISPGFPLOT is true, two pgfplot-files are created with the names
%  pNAME.pgfplot and TNAME.pgfplot, respectively. With ISQPLOT being true, a q-z
%  plot is created. VPOS gives the vertical position of the figure's bottom on
%  the screen:

%fl.info.T1, .p1, .p2, .substance, .membrane, .fmodel, .flsetup
%fl.sol.T2, .T1, .len, .colors
%fl.calc.psat1, .pK1, .n

if nargin < 5
  vpos = 535;
  if nargin < 4
    isqplot = false;
    if nargin < 3
      ispgfplot = false;
      if nargin < 2
	name = [];
      end
    end
  end
end


% The number of line-parts.
nflow = -fl.sol.len;
L = fl.info.membrane.L;
isfront = strcmp('13',fl.sol.states(1:2));
if isfront
  % scale the front boundary layer, it may not always be present
  z3 = fl.flow(nflow).z(1);
  zscale = fl.sol.zscale;
  % Aus flow12.m
  %    z3 = FL.flow(-FL.sol.len).z(1),  zscale = FL.sol.zscale,
  %    z = z3 + zscale * log( (Fl.flow(-FL.sol.len).z-z3)/zscale + 1 ).
  %z = [fl.flow(1:end-1).z z3+zscale*log((fl.flow(nflow).z-z3)/zscale+1)]/L;
  % the last point might become, non-dimensionalized, less than zero; correct
  % this to (minus) infinity
  if (-fl.flow(nflow).z(end)+z3) < zscale
    zfront = z3 + zscale*log((fl.flow(nflow).z-z3)/zscale+1);
  else
    zfront = z3 + zscale*log((fl.flow(nflow).z(1:end-1)-z3)/zscale+1);
    zfront = [zfront 10*zfront(end)];
  end
else
  zfront = fl.flow(nflow).z;
  zscale = 0;
end

% Flattened distributions.
% only used for the plot limits, immediately below
Tflat = [fl.flow(1:end).T];
pflat = [fl.flow(1:end).p];
if isqplot, qflat = [fl.flow(1:end).q]; end

% The temperature range
margin = 2; step = 2;
Tmax = ceil( (max(Tflat)+margin)/step ) * step;
Tmin = floor( (min(Tflat)-margin)/step ) * step;

% The pressure range [Pa]
margin = 1e4; step = 2e4;
pmax = ceil( (max(pflat)+margin)/step ) * step;
pmin = floor( (min(pflat)-margin)/step ) * step;
pmax = pmax / 1e5;
pmin = pmin / 1e5;

% q-range
if isqplot
  margin = 5; step = 10;
  qmax = ceil( (max(qflat)+margin)/step ) * step;
  qmin = floor( (min(qflat)-margin)/step ) * step;
end

% z-range
step = 0.2; margin = 0.1;
Tzmin = floor((-3*zscale/L-margin)/step) * step;
pzmin = -0.1;
zmax = 1.1;

%%%	THE PLOTS	%%%
Tname = [name 'Tz'];
pname = [name 'pz'];
if isqplot, qname = [name 'qz']; end

% pgfplots
if ispgfplot
  % T-z diagram
  rangestr = sprintf(' xmin =%.1f, xmax = %.1f, ymin = %.0f, ymax = %.0f,\n',...
	Tzmin, zmax, Tmin, Tmax);
  Tid = beginpgfplot(Tname,['xlabel = {$z/L$}, ylabel = {$T$ [K]},\n' rangestr ...
    ' legend style={at={(0.97,0.07)},anchor=south east,cells={anchor=west}},\n'...
    ' y label style = {rotate=-90,xshift=-10bp}, width=8cm, height=6cm']);

  addcoords(Tid,zfront'/L,fl.flow(nflow).T',...
	'mark=*,mark options={scale=0.6},solid,thick');
  for i = nflow-1:-1:1
    addcoords(Tid,fl.flow(i).z'/L,fl.flow(i).T',...
	'mark=*,mark options={scale=0.6},solid,thick');
  end
  endpgfplot(Tid);

  % p-z diagram
  rangestr = sprintf(' xmin =%.1f, xmax = %.1f, ymin = %.1f, ymax = %.1f,\n',...
	pzmin, zmax, pmin, pmax);
  pid = beginpgfplot(pname, [...
    'xlabel = {$z/L$}, ylabel = {$p$ [bar]},\n' rangestr ...
    ' legend style={at={(0.97,0.07)},anchor=south east,cells={anchor=west}},\n'...
    ' y label style = {rotate=-90,xshift=-10bp}, width=8cm, height=6cm']);

  addcoords(pid,zfront'/L,fl.flow(nflow).p'/1e5,...
	'mark=*,mark options={scale=0.6},solid,thick');
  for i = nflow-1:-1:1
    addcoords(pid,fl.flow(i).z'/L,fl.flow(i).p'/1e5,...
	'mark=*,mark options={scale=0.6},solid,thick');
  end
  endpgfplot(pid);
else
% matlab-plots
  % T-z diagram
  figure('Name',['T-z diagram:  ' Tname],'OuterPosition',[100 vpos 562 504]);
  %plot(zflat,T,'k*');
  plot(zfront/L,fl.flow(nflow).T,'Color',fl.flow(nflow).color,...
	'LineStyle','-','Marker','*');
  xlim([Tzmin zmax]);
  xlabel('z/L');
  ylabel('T [K]');
  for i = nflow-1:-1:1
    line(fl.flow(i).z/L,fl.flow(i).T,'Color',fl.flow(i).color,...
	'Marker','*','LineStyle','-');
  end
  % p-z diagram
  figure('Name',['p-z diagram:  ' pname],'OuterPosition',[720 vpos 562 504]);
  plot(zfront/L,fl.flow(nflow).p/1e5,'Color',fl.flow(nflow).color,...
	'LineStyle','-','Marker','*');
  xlim([pzmin zmax]);
  ylim([pmin pmax]);
  xlabel('z/L');
  ylabel('p [bar]');
  % hold on; plot..; hold off ginge auch, ohne xlim, ylim setzen zu müssen
  for i = nflow-1:-1:1
    line(fl.flow(i).z/L,fl.flow(i).p/1e5,'Color',fl.flow(i).color,...
	'Marker','*','LineStyle','-');
  end
  % q-z diagram
  if isqplot
    figure('Name',['q-z diagram:  ' qname],'OuterPosition',[1340 vpos 562 504]);
    plot(zfront/L,fl.flow(nflow).q,'Color',fl.flow(nflow).color,...
	  'LineStyle','-','Marker','*');
    xlim([Tzmin zmax]);
    ylim([qmin qmax]);
    xlabel('z/L');
    ylabel('q [W/m2]');
    % hold on; plot..; hold off ginge auch, ohne xlim, ylim setzen zu müssen
    for i = nflow-1:-1:1
      line(fl.flow(i).z/L,fl.flow(i).q,'Color',fl.flow(i).color,...
	  'Marker','*','LineStyle','-');
    end
  end
end
end
