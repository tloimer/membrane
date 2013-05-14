function Tsplot(fl,name,ispgfplot)
%TSPLOT     Plot a T-s Diagram.
%  TSPLOT(NAME,FLOWSTRUCT,ISPGFPLOT) plots a T-s diagram. If ISPGFPLOT is true,
%  a pgfplot-file with the name NAME.pgfplot is created.

s = fl.info.substance;
% the reference point is the lowest on the liquid line, (Tr,ps(Tr),0).
Tr = fl.sol.T2;
T1 = fl.info.T1;
T2 = fl.sol.T2;
p1 = fl.info.p1;
p2 = fl.info.p2;

% the region of interest
dT = T1 - T2;
Tminr = T2 - dT/2;
Tmaxr = T1 + dT/2;
Trange = Tminr : dT/10 : Tmaxr+dT/100; % to prevent rounding errors
nTrange = size(Trange,2);

% the lines of saturated vapor and liquid
fs = fl.info.flsetup;
psrange = Trange;
pKrange = Trange; pKcaprange = Trange;
for i = 1:nTrange
  psrange(i) = s.ps(Trange(i));
  pKrange(i) = psrange(i)*fs.pkps(Trange(i));
  pKcaprange(i) = pKrange(i) - fs.curv*s.sigma(Trange(i));
end

sgrange = s.s(Trange,psrange,1,Tr);
slrange = s.s(Trange,psrange,0,Tr);
sKrange = s.s(Trange,pKrange,1,Tr);
scrange = s.s(Trange,pKcaprange,0,Tr);

% bounds of plotted area
% sgrange(1) is the entropy of vaporization at Tminr
smin = floor(-sgrange(1)/3/200)*200;
smax = ceil(1.5*sgrange(1)/200)*200;
Tmin = floor(Tminr);
Tmax = ceil(Tmaxr);

% the isobars p1 and p2
Ts1 = s.Ts(p1);
Ts2 = s.Ts(p2);
if Ts1 >= Tminr,
  Tl1range = [Tminr:dT/10:Ts1 Ts1];
else
  Tl1range=[];
end
Tv1range = [max(Ts1,Tminr):dT/10:Tmaxr];
if Ts2 >= Tminr,
  Tl2range = [Tminr:dT/10:Ts2 Ts2];
else
  Tl2range=[];
end
Tv2range = [max(Ts2,Tminr):dT/10:fl.info.T1];
sp1range = [s.s(Tl1range,p1,0,Tr) s.s(Tv1range,p1,1,Tr)];
sp2range = [s.s(Tl2range,p2,0,Tr) s.s(Tv2range,p2,1,Tr)];

% and the isenthalpic line from T1, p1 to max(p1-3*p12,100 Pa)
p12 = p1 - p2;
prange = p1:-p12/12:max(100,p1-3*p12);
prange = prange(2:end);
Thrange = [T1 s.intjt(T1,p1,prange)];
shrange = s.s(Thrange,[p1 prange],1,Tr);
if s.ps(T1) > p1
  % prange contains decreasing pressures
  prange = min(s.ps(T1),p1+p12):-p12/12:p1;
  % intjt must integrate from p1 away; turn around prange
  Thuprange = s.intjt(T1,p1,prange(end:-1:1));
  % temperatures in Thrange are now increasing
  % turn around Thrange
  Thuprange = Thuprange(end:-1:1);
  shrange = [s.s(Thuprange,prange,1,Tr) shrange];
  Thrange = [Thuprange Thrange];
end

% Calculate the entropy for the path
nflow = abs(fl.sol.len);
Tpath = cell(nflow);
spath = cell(nflow);
for i = 1:nflow
  Tpath{i} = fl.flow(i).T;
  if all(fl.flow(i).a == 1) || all(fl.flow(i).a == 0)
    spath{i} = s.s(Tpath{i},fl.flow(i).p,fl.flow(i).a,Tr);
  else
    % pre-allocate memory
    spath{i} = Tpath{i};
    for ii = 1:length(fl.flow(i).T)
      Tii = fl.flow(i).T(ii);
      pK = fs.pkps(Tii)*s.ps(Tii);
      pKcap = pK - fs.curv*s.sigma(Tii);
      % ouh, must calculate the mass fraction
      x = fl.info.fmodel.x(fl.flow(i).a(ii),s.v(Tii,pK),1/s.rho(Tii));
      spath{i}(ii) = (1-x)*s.s(Tii,pKcap,0,Tr) + x*s.s(Tii,pK,1,Tr);
    end
  end
end

% Start plotting
if ispgfplot
  rangestr = sprintf(' xmin = %.0f, xmax = %.0f, ymin = %.0f, ymax = %.0f,\n',...
	smin, smax, Tmin, Tmax);
  sid = beginpgfplot(name,['xlabel={$s-s_0$ [J/kg\\,K]}, ylabel={$T$ [K]},\n'...
    rangestr ...
    ' legend style={at={(0.97,0.07)}, anchor=south east, cells={anchor=west},\n'...
    '\tdraw = white},\n'...
    ' y label style = {rotate=-90,xshift=1bp}, width=8cm, height=6cm']);

  % Saturated liquid
  fprintf(sid,'%% Saturated liquid\n');
  addcoords(sid,slrange',Trange','blue!30!white, solid, thick');
  % Saturated vapor
  fprintf(sid,'%% Saturated vapor\n');
  addcoords(sid,sgrange',Trange','orange!90!red, solid, thick');
  % pk and pk-pcap dotted lines % -- with pgfplots-injection!
  fprintf(sid,'%% pK - pcap\n');
  addcoords(sid,scrange',Trange','thin, dotted] plot[forget plot');
  fprintf(sid,'%% pK\n');
  addcoords(sid,sKrange',Trange','thin, dotted');
  % isobars p1 and p2
  fprintf(sid,'%% isobar p = p1\n');
  addcoords(sid,sp1range',[Tl1range Tv1range]','black!70!white');
  fprintf(sid,'%% isobar p = p2\n');
  addcoords(sid,sp2range',[Tl2range Tv2range]','black!70!white');
  % isenthalpic line
  fprintf(sid,'%% isenthalpic line, h = h1\n');
  addcoords(sid,shrange',Thrange','black!70!white, thin, dashed');

  % the path
  fprintf(sid,'%% path of the process\n');
  for i = 1:nflow
    addcoords(sid,spath{i}.',Tpath{i}.',' very thick, solid');
    %addcoords(sid,spath{i}(1),Tpath{i}(1),'mark=*, none');
    %addcoords(sid,spath{i}(end),Tpath{i}(end),'mark=o, none');
  end
  endpgfplot(sid,['saturated liquid, saturated vapor, curved interface, ' ...
	'$p=p_1$, $p=p_2$, $h=h_1$, path of the process']);
else % matlab-plots
  Tsh = figure('Name',['T-s diagram: ' name]); %,'NextPlot','add');
  axes('NextPlot','add');
  % Saturated liquid
  plot(slrange,Trange,'Linewidth',2,'Color',[0.5 .5 1]);
  xlim([smin smax]);
  ylim([Tmin Tmax]);
  xlabel('s-s_0 [J/kgK]');
  ylabel('T [K]');
  % Saturated vapor
  line(sgrange,Trange,'Linewidth',2,'Color',[1 0.6 0.2]);
  % pk and pk-pcap dotted lines
  plot(scrange,Trange,'k:',sKrange,Trange,'k:')
  % isobars p1 and p2
  line(sp1range,[Tl1range Tv1range],'Color',[0.6 0.6 0.6]);
  line(sp2range,[Tl2range Tv2range],'Color',[0.6 0.6 0.6]);
  % isenthalpic line
  line(shrange,Thrange,'Color','k','LineStyle','--');

  for i = 1:nflow
    line(spath{i},Tpath{i},'Color','k','LineWidth',2);
    line(spath{i}(1),Tpath{i}(1),'Color','k','Marker','o','MarkerSize',6);
    line(spath{i}(end),Tpath{i}(end),'Color','k','Marker','.','MarkerSize',18);
  end
  set(Tsh,'NextPlot','new');
end
