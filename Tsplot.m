function Tsplot(fs)
%TSPLOT     Plot a T-s Diagram.
%  TSPLOT(FLOWSTRUCT)
%
%  Calls INTG, INT2PH, INTP, TEMP, Q_M, X, R, PS, TS, CPL, K.
%  Subfunctions: S, SG, SL, H, HG, HL, ISOP, ISOH.
%  See also FLOWSTRUCT.

% the reference point is the lowest on the liquid line, (Tr,ps(Tr),0).
Tr = fs.sol.Te;

% the region of interest
dT = fs.info.T0 - fs.sol.Te;
Tmin = fs.sol.Te - dT/2;
Tmax = fs.info.T0 + dT/2;
Trange = Tmin : dT/15 : Tmax;

% the complete path; add the initial point of saturated vapor
T = [fs.flow(:).T fs.info.T0];
p = [fs.flow(:).p fs.info.p0];
xx = [fs.flow(:).x 1];

% corner points
len = abs(fs.sol.len);
% memory allocation + last point assignment
Tcp(2*len+1) = fs.info.T0;
scp(2*len+1) = s(fs.info.T0,fs.info.p0,1,Tr);
for i = 1:len
  n = 2*i;
  Tcp(n-1) = fs.flow(i).T(1);
  Tcp(n) = fs.flow(i).T(end);
  scp(n-1) = s(Tcp(n-1),fs.flow(i).p(1),fs.flow(i).x(1),Tr);
  scp(n) = s(Tcp(n),fs.flow(i).p(end),fs.flow(i).x(end),Tr);
end

% the first point at the front (after condensation)
Tf = T(end-1);
pf = p(end-1);
xf = xx(end-1);

% and calculate the entropies
sy = s(T,p,xx,Tr);

% additional stuff in the diagram
prange = ps(Trange);

% the liquid and vapor line
liq = s(Trange,prange,0,Tr);
vap = s(Trange,prange,1,Tr);

% and some mass fraction lines
%x01 = s(Trange,prange,0.1,Tr);
x02 = s(Trange,prange,0.2,Tr);
x05 = s(Trange,prange,0.5,Tr);
x08 = s(Trange,prange,0.8,Tr);

% isobars at p0 and pe (fs.sol.pe, not fs.info.p0-fs.info.dp)
[Tp0 sp0] = isop(Trange,fs.sol.p0,Tr);
[Tpe spe] = isop(Trange,fs.sol.pe,Tr);
[Tv sv] = isop(Trange,fs.flow(1).p(end),Tr);

% isenthalps at (T0,p0) and at the front of the membrane
[Th0 sh0] = isoh(Trange,fs.info.T0,fs.info.p0,1,Tr);
[Thf shf] = isoh(Trange,Tf,pf,xf,Tr);

% first plot all the underlying stuff
plot(liq,Trange,vap,Trange,'LineWidth',2,'Color',[0.8 0.8 0.8]);
line(0,Tr,'Marker','o','Color','k');
line([x02;x05;x08],Trange,'Color','k','LineStyle','-.');
line([sp0;sv;spe]',[Tp0;Tv;Tpe]','Color','k','LineStyle','--');
line(sh0,Th0,'Color','k','LineStyle','--');
line(shf,Thf,'Color','k','LineStyle','--');
%line(sv,Tv,'Color','r','LineStyle','--');
%line([sp0 spe],[Tp0 Tpe],'Color','r','LineStyle','--');

% some points are labeled
%text(sy(end),fs.info.T0,'1','VerticalAlignment','bottom');
%text(sy(1),T(1),'2','VerticalAlignment','bottom');

% add the points along the linear path
Tl(6) = 0; sl(6) = 0;
Tl(1) = fs.info.T0; sl(1) = s(fs.info.T0,fs.info.p0,1,Tr);
Tl(2) = fs.info.T0; sl(2) = s(fs.info.T0,fs.info.p0,fs.lin.x3,Tr);
% in case of 2ph-flow, the following point is identical to Tl(2)
if ~isempty(fs.lin.T4)
  Tl(3) = fs.lin.T4;  sl(3) = s(Tl(3),fs.info.p0,fs.lin.x3,Tr);
else
  %Tl(3) = []; sl(3) = [];
  Tl(3) = Tl(2); sl(3) = sl(2);
end
Tl(4) = fs.lin.Te;  sl(4) = s(Tl(4),ps(Tl(4)),fs.lin.x3,Tr);
Tl(5) = Tl(4);      sl(5) = s(Tl(4),ps(Tl(4)),1,Tr);
Tl(6) = Tl(4);      sl(6) = s(Tl(4),fs.info.p0-fs.info.dp,1,Tr);
line(sl,Tl,'Color','y','LineStyle','-','LineWidth',0.5);%,'Marker','x');

% at last the nonlinear path
line(sy,T,'Color','k');
line(scp,Tcp,'Color','k','LineStyle','none','Marker','x');

function sy = s(T,p,x,Tr)
%S          Specific entropy [J/kgK].
% S(T,P,X,TR) returns the difference of the entropies S(T,P,X) -
% S(TR,PS(TR),0). The reference state lies on the liquid line.

% return empty for empty input
%if isempty(T) | isempty(p) | isempty(x)
  %error('s called with empty T,p,x');
  %sy = [];
  %return
%end

%--------- delta s is calculated along an isobar
% input check
% we have to go step by step

% make all vectors
lT = length(T);
lp = length(p);
lx = length(x);
i = max([lT lp lx]);
if i>1
  if lT == 1
    T(1:i) = T;
  end
  if lp == 1
    p(1:i) = p;
  end
  if lx == 1
    x(1:i) = x;
  end
end

sy(1:i) = 0;
prec = 1 + 1e-8;
for j = 1:i
  switch x(j)
    case 1
      if p(j)>prec*ps(T(j))
        error(sprintf('State of subcooled vapor at index %d',j));
      end
      % vapor contribution
      T1 = Ts(p(j));
      sy(j) = sg(T1,T(j),p(j)) + r(T1)/T1;
      T(j) = T1;
  	% make sure p=ps(T1)
  	p(j) = ps(T1);
    case 0
      if prec*p(j)<ps(T(j))
        error('State of superheated liquid');
      end
    otherwise
      if x(j)>1 | x(j)<0
        error(sprintf('Nonsense input: x = %.3g.',x(j)));
      else % 0<x<1
        if p(j)~=ps(T(j))
          error('Something strange happened: p not equal to ps(T)');
        end
        sy(j) = x(j)*r(T(j))/T(j);
      end
  end
end

% now the state is in the liquid, or on the evap. line
% liquid contribution
sy = sy + sl(Tr,T,ps(Tr),p);

function syg = sg(T1,T2,p1,p2)
%SG         Difference of the specific entropy of the vapor [J/kgK].
% SG(T1,T2,P1,P2) returns the entropy difference S2(T2,P2) - S1(T1,P1),
% based on the thermal equation of state in the vapor region.
%
% SG(T1,T2,P) returns the entropy difference S2(T2,P) - S1(T1,P).
%
% Calls MOLM, TDBDT.

[R M] = molm;
[B TdBdT] = TdbdT(T1);
dbdT1 = TdBdT./T1;
[B TdBdT] = TdbdT(T2);
dbdT2 = TdBdT./T2;

lT = length(T1);
dsT = zeros(1,lT);
dcp = zeros(1,lT);

for i = 1:lT
  if T2(i)~=T1(i)
    dcp(i) = cpid(T2(i)) - cpid(T1(i));
    dsT(i) = (cpid(T1(i))-dcp.*T1(i)./(T2(i)-T1(i))).*log(T2(i)/T1(i)) + dcp(i);
  end
end

if nargin == 3
  % isobaric
  syg = dsT - p1.*(dbdT2-dbdT1)/(1000*M);
elseif nargin == 4
  syg = dsT + ( (dbdT1.*p1-dbdT2.*p2)/1000 - R*log(p2./p1) )/M;
end

function syl = sl(T1,T2,p1,p2)
%SL         Difference of specific entropy for the liquid [J/kgK].
% SL(T1,T2,P1,P2) returns the entropy difference S2(T2,P2) - S1(T1,P1)
% in the liquid region. An incompressible liquid is assumed, but the
% density RHO depends on the temperature T.
%
% Calls MOLM, TDBDT.

syl = intcpl_T(T1,T2) - (p2-p1)./(rho(T1).*T1);

function hy = h(T,p,x,Tr)
%H          Specific enthalpy [J/kg].
% H(T,P,X,TR) returns the difference of the enthalpy H(T,P,X) -
% H(TR,PS(TR),0). The reference state lies on the liquid line.

% input check
% we have to go step by step

% make all input arguments vectors
lT = length(T);
lp = length(p);
lx = length(x);
i = max([lT lp lx]);
if i>1
  if lT == 1
    T(1:i) = T;
  end
  if lp == 1
    p(1:i) = p;
  end
  if lx == 1
    x(1:i) = x;
  end
end

hy(1:i) = 0;
prec = 1 + 1e-8;
for j = 1:i
  switch x(j)
    case 1
      if p(j)>prec*ps(T(j))
        error(sprintf('State of subcooled vapor at index %d',j));
      end
      % vapor contribution
      T1 = Ts(p(j));
      hy(j) = hg(T1,T(j),p(j)) + r(T1);
      T(j) = T1;
  	% make sure p=ps(T1)
  	p(j) = ps(T1);
    case 0
      if prec*p(j)<ps(T(j))
        error('State of superheated liquid');
      end
    otherwise
      if x(j)>1 | x(j)<0
        error(sprintf('Nonsense input: x = %.3g.',x(j)));
      else % 0<x<1
        if p(j)~=ps(T(j))
          error('Something strange happened: p not equal to ps(T)');
        end
        hy(j) = x(j)*r(T(j));
      end
  end
end

% now the state is in the liquid, or on the evap. line
% liquid contribution
hy = hy + hl(Tr,T);

function hyg = hg(T1,T2,p1,p2)
%HG         Difference of the specific enthalpy of the vapor [J/kg].
% HG(T1,T2,P1,P2) returns the entropy difference H2(T2,P2) - H1(T1,P1),
% based on the thermal equation of state in the vapor region.
%
% HG(T1,T2,P) returns the entropy difference H2(T2,P) - H1(T1,P).
%
% Calls CPID, DHDP.

if nargin==3
  p2=p1;
end
if length(T1)~=1 | length(T2)~=1 | length(p1)~=1 | length(p2)~=1
  error ('not vectorized.')
end

% integrate cp dT
if T2==T1
  hyg = (p2-p1)*dhdp(T1);
else
  dcp = ( cpid(T2) - cpid(T1) )/(T2-T1);
  hyg = (cpid(T1)-dcp*T1)*(T2-T1) + dcp*(T2^2-T1^2)/2 ...
    + p2*dhdp(T2) - p1*dhdp(T1);
end
%hyg = hyg + p2*dhdp(T2) - p1*dhdp(T1);

function hyl = hl(T1,T2)
%HL         Difference of specific enthalpy for the liquid [J/kg].
% HL(T1,T2) returns the enthalpy difference H2(T2) - H1(T1)
% in the liquid region.
%
% Calls INTCPL.

hyl = intcpl(T1,T2);

function [Tp,sp] = isop(Trange,p,Tr)
%ISOP       Construct an isobar.
% [TP,SP] = ISOP(TRANGE,P) returns vor a vector TRANGE a vector TP and a
% vector SP with the corresponding entropies.

T = Ts(p);
if T>Trange(1) & T<Trange(end)
  Tliq = Trange(find(ps(Trange)<=p)); % in case one ps(Trange)==ps, Tp
  Tvap = Trange(find(ps(Trange)>=p)); % length(Tp) is length(Trange)+1
  lliq = length(Tliq);
  lvap = length(Tvap);
  Tp = [Tliq(1:lliq-1) T T Tvap(2:lvap)];
  x(1:lliq) = 0;
  x(lliq+1:lliq+lvap) = 1;
else
  Tp = Trange;
  if T<Trange(1)
    x = 1;
  else
    x = 0;
  end
end

sp = s(Tp,p,x,Tr);

function [Th,sh] = isoh(Trange,T0,p0,x0,Tr)
%ISOH       Construct an isenthalpe.

% calculate with reference to (T0,ps(T0))
hy = h(T0,p0,x0,T0);

% is the isenthalp in the graph region?
if hy < hl(Tr,Trange(1))
  Th = [];
  sh = [];
  warning('this isenthalp is not in the diagram region');
  return
end

% care for the part in the liquid region
if hy < hl(T0,Trange(end))
  Tend = fzero(@hldif,[Trange(1) Trange(end)],optimset('fzero'),hy,T0);
  Th = Trange(find(Trange<Tend));
  lT = length(Th);
  % and add two points for a horizontal line
  Th(lT+1) = Tend;
  Th(lT+2) = Tend;
  sh(lT+1) = s(Tend,ps(Tend),0,Tr);
  sh(lT+2) = 0;
else
  lT = length(Trange);
  Th = Trange;
  sh = zeros(1,lT);
end

ilast = 0; % to delete points
for i=1:lT
  hvap = r(T0)+hg(T0,Th(i),ps(T0),ps(Th(i)));
  hliq = hvap - r(Th(i));
  %hvap = hliq + r(Th(i));
  if hvap>=hy
    x = (hy-hliq)/(hvap-hliq);
    sh(i) = s(Th(i),ps(Th(i)),x,Tr);
  else
    x=1;
    if hdif(0,T0,Th(i),hy)<0
    %hdif(0,T0,Th(i),hy)>hdif(ps(Th(i)),T0,Th(i),hy)
    % the isenthalp becomes horizontal
      ilast = i;
    else
      p = fzero(@hdif,[0 ps(Th(i))],optimset('fzero'),T0,Th(i),hy);
      sh(i) = s(Th(i),p,1,Tr);
    end
  end
end
Th(1:ilast) = [];
sh(1:ilast) = [];

function dh = hdif(p,T0,T,hy)
dh = r(T0) + hg(T0,T,ps(T0),p) - hy;

function dh = hldif(T,hy,Tr)
dh = hy - hl(Tr,T);
