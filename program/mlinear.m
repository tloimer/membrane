function [m,fl] = mlinear(p1,p2,T1,theta,s,mem,f)
%MLINEAR    Mass flux from linear theory through a homogeneous membrane.
%  MLINEAR(P1,P2,T1,THETA,S,M,F) calculates the mass flux [kg/m2s] of
%  substance S through membrane M according to a linear theory. The contact
%  angle is given by THETA and F is the two-phase flow model. For P1 =
%  PSAT(T1), the mass flux is calculated as given in JMS07 [Loimer,
%  J. Membr. Sci. 301, pp. 107-117, 2007]. For wetted systems and
%  PK - N*(P1-P2) < P1 < PSAT, the mass flux is obtained from an linear
%  interpolation between MGAS and MLINEAR(P1=PSAT), see JMS11 [Loimer, J.
%  Membr, Sci. 383, pp. 104-115, 2011]. P1 [Pa] is the upstream pressure,
%  P2 [Pa] the downstream pressure, and T1 [K] is the upstream temperature.
%  The substance S, membrane M and two-phase flow model F are provided via
%  structs, see SUBSTANCE, MEMBRANE and FMODEL.
%
%  [M,FL] = MLINEAR(P1,P2,T1,THETA,S,M,F) writes the solution for P1 =
%  PSAT(T1), ignoring P1, to the FLOWSTRUCT FL. Use PTPLOT, PTZPLOTS and
%  TSPLOT to plot the solution stored in FL.
%
%  See also FMODEL, MEMBRANE, SUBSTANCE, PTPLOT, PTZPLOTS, TSPLOT and
%  MLINEAR>MLINPSAT.

% Some input sanitizing.
if s.ps(T1) < p1
  error([upper(mfilename)...
	': The upstream state is a liquid. This is not implemented.']);
elseif p2 < 0.
  error([upper(mfilename)...
	': The downstream pressure is negative. That is not possible.']);
elseif T1 < 0.
  error([upper(mfilename)...
	': The upstream temperature is below absolute zero. Impossible.']);
end

% THETA is in degree - calculate costheta.
costheta = cos(theta*pi/180);
p12 = p1 - p2;

% Necessary variables.
[ps dps] = s.ps(T1);
R = s.R;  vliq = 1/s.rho(T1); vgas = s.v(T1,(p1+p2)/2);
mugas = s.mug(T1); sigma = s.sigma(T1); mujt = s.jt(T1,p1);

dia=mem.dia; beta = mem.beta; L = mem.L;
kap = mem.kappa; curv = mem.fcurv(costheta);

% Calculate the apparent kinematic viscosity of the gas.
facKn = 3*sqrt(pi/(8*R*T1));
nubar = mugas*vgas;
nuapp=1/(1/nubar + beta*facKn/dia);

% Calculate mgas, mgas = kappa*(p1-p2)/(nuapp*L).
if abs(costheta)>=1e-3
  % eq. (6)
  pcap = curv*sigma;
  % eq. (7)
  pk_ps = exp(-pcap*vliq/(R*T1));
  pk = pk_ps*ps;
  % eq. (11), truncated at the first term; eq. (13)
  n = mujt*pk_ps*dps;
else % abs(costheta)<=1e-3
  costheta = 0;
  pcap = 0; pk = ps; n = mujt*dps;
end
mgas = kap*p12/(nuapp*L);

pref = pk - n*p12;

if p1 <= pref
  m = mgas;
else
  fl = mlinpsat(T1,p12,theta,s,mem,f);
  % m = ( mgas*(ps-p1) + fl.lin.m*(p1-pref) ) / (ps-pref);
  % probably more accurate, see mnum.m, line 121
  m = mgas + (fl.lin.m-mgas)*(p1-pref)/(ps-pref);
end
%%% END ASYM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END ASYM %%%

%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

function fl = mlinpsat(T1,p12,theta,s,mem,f) %------------------------- mlinpsat
%MLINPSAT   Mass flux from linear theory for P1 = PSAT(T1) [kg/m2s].
%  MLINPSAT(T1,P12,THETA,S,M,F) calculates the mass flux according to
%  linear theory [Loimer, J. Membr. Sci. 301, pp. 107-117, 2007]. The
%  upstream pressure is PSAT(T1), with T1 given in Kelvin. A pressure
%  difference P12 [Pa] is applied. The substance S, membrane M and
%  two-phase flow model F are provided via structs, see SUBSTANCE, MEMBRANE
%  and FMODEL. The output is a (obsolete) flowstruct FL. The mass flux is
%  reported in FL.LIN.M (or FL.INFO.M) [kg/m2s].
%
%  See also FLOWSTRUCT, SUBSTANCE, MEMBRANE, FMODEL,
%           MSTACKSTRUCT>SINGLEMSTOFL.
%
%  Adapted from MLINEAR.M in, e.g., 12potsdam, 11jms/matlab3.

%% Input processing.
% THETA is in degree - calculate costheta.
costheta = cos(theta*pi/180);

%% Setup of material properties.
% Calls to SUBSTANCE.M and TOPOLOGY.M. The 2ph-model as well as membrane
% material properties (e.g., km, epsilon, beta) are set here.

% properties, fluid
[ps dps] = s.ps(T1); p1 = ps;
R = s.R;  vliq = 1/s.rho(T1); vgas = s.v(T1,p1); hvap = s.hvap(T1);
cpgas = s.cpg(T1,p1); muliq = s.mul(T1); mugas = s.mug(T1);
kgas = s.kg(T1); kliq = s.kl(T1); sigma = s.sigma(T1); mujt = s.jt(T1,p1);
% cpliq = s.cpl(T1);
% auxiliary variables, depending on substance only
nuliq = muliq*vliq; % nugas=mugas*vgas;
% In linear theory, the enthalpy of vaporization cannot depend on curvature.
rk = hvap;

% properties, membrane (material and structural (=topology) properties)
km = mem.km; dia=mem.dia; epsilon = mem.epsilon; beta = mem.beta;
L = mem.L; kap = mem.kappa; curv = mem.fcurv(costheta);

% 2ph-flow model
% Functions to calculate 2ph-flow properties. Only homogeneous flow ( = 'plug')
% is used, evaluated at the upstream state T1, p1, hence the functions are given
% here.
fxdot = @(a) f.xdot(a,vgas,vliq);
fx = @(a) f.x(a,vgas,vliq);
fnu2ph = @(a) f.nu2ph(a,vgas,vliq,mugas,muliq);
fk2ph = @(a) f.k2ph(a,epsilon,km,kgas,kliq);
%fmodel=f.name; % for struct fl needed

% variables depending on more than one of the called property functions
% (SUBSTANCE, TOPOLOGY, 2ph-model).
kml = fk2ph(0);

%% Calculation of dependent variables - pk, pcap, etc.
p2 = p1 - p12;
% eq. (8)
T12 = mujt*p12; T2 = T1 - T12;

% TL. p. 109, above eq. (4)
% From Bird, Stewart, Lightfoot, p.20: mean free path
% lambda = 3 nu sqrt(pi/(8RT)),
% hence lambda = facKn*nu and Kn = facKn*nu/dia
facKn = 3*sqrt(pi/(8*R*T2));

% if costheta~=0
if abs(costheta)>=1e-3
  % eq. (6)
  pcap = curv*sigma;
  % eq. (7)
  pk_ps = exp(-pcap*vliq/(R*T1));
  pk = pk_ps*ps;
  p1pk = p1-pk;
  % eq. (11), truncated at the first term
  dpk = pk_ps*dps;
  % eq. (13)
  n = mujt*dpk;
  if n>=1
    error([upper(mfilename) ': state 2 is a 2ph-mixture: n >= 1, n = %f'],n);
  end

  if costheta >= 0
    % eq. (12)
    Cc = p12*(1-n)/p1pk; % Cc > 0; Cc = Ccc;
  else
    % eq. (14)
    Cc = n*p12/pcap; % Cc < 0; Cc = Ccap;
  end
else % abs(costheta)<=1e-3
  %disp(sprintf('theta corrected to 90°, cos(theta) = %g',costheta));
  costheta = 0; Cc = [];
  pcap = 0; pk = ps; p1pk = 0; dpk = dps; n = mujt*dps;
end

% eq. (10): kapc is really kappaK(kappa). Eq. (10) is implicit, because dp_K/dT
% is a function of kappa. Hence, one would have to distinguish between the
% function kappaK(kappa) and an unique value of kappaK, which is given by
% kappaK(kappaK). For calculation of the mass fluxes, the function kappaK(kappa)
% must be used.
kapc = nuliq*kml/(dpk*hvap);
kkc=kap/kapc;

%% Calculate the apparent vapor viscosity.
% The kinematic viscosity of the vapor is calculated assuming ideal-gas law for
% the vapor and applying a correction for molecular flow.

% eq. (45)
p9=pk-n*p12;
% nubar, see TL, p. 110, top three lines
% $$ \bar\nu = \mu (v_9 + v_2) / 2. $$
% with Boyle-Marriott, p1 v1 = pmean vmean, pmean = (p2+p9)/2, T = T1
nubar = mugas*vgas*2*p1/(p2+p9);
% DEBUG
% probably, take mugas and vgas at T2:
%mug2 = s.mug(T2); vg2 = s.v(T2,p2);
% nubar = mug2*vg2*2*p2/(p2+p9);
%tmp = mug2*vg2*2*p2/(p2+p9);
%disp(sprintf('nubar with downstream-values, nubar(T2,p2)/nubar(T1,p2) = %f',...
% tmp/nubar));
% END DEBUG
% eq. (3), Kn = nubar*facKn/dia
nuvap=1/(1/nubar + beta*facKn/dia);

% p9 stimmt fuer non-wetting nicht ganz - aber gut genug
% Route zu w3, w1, w2 checked.
%disp(sprintf('kkc = %g,  b = %g nm,  n = %g,  Cc = %g,  kml = %g',...
%	kkc,b*1e9,n,Cc,kml));
%disp(sprintf('p1pk = %g Pa,  p9 = %g bar,  Kn = %g',p1pk,p9/1e5,...
%  nubar*facKn/dia));
%disp(sprintf('nubar = %g uPas,  nuvap = %g uPas',nubar*1e6,nuvap*1e6));

%% Prepare output.
% Construct the output struct, see also FLOWSTRUCT.

% lin.Te is the temperature at the evaporation front
fl = struct('info',struct('kap',kap,'m',[],'C',Cc,'kapc',kapc,'kapf',[],...
  'L',L,'T0',T1,'T1',T1,'p0',p1,'p1',p1,'dp',p12,'p2',p2,'ph',f.name,...
  'flsetup',flowsetup(T2,T2+1.2*T12,theta,s,mem,f)),...
  'sol',struct('len',[],'a3',[],'q3',[],'T0',T1,'Te',T2,'T2',T2,'p0',p1,...
    'pe',p2,'de',[],'df',[],'states','--'),...
  'flow',struct('z',{},'T',{},'p',{},'q',{},'a',{},'x',{},'color',{}),...
  'lin',struct('m',[],'Te',T2,'T4',[],'deL',[],'dfL',[],'a3',[],'x3',[]));

%% Determine the type of flow according to the flow map, Fig. 1.

% costheta==0:  Schneider's solutions_____
%   kkc<=1		s1	liquid film - liq. flow - vapor or 2ph
%   kkc>1		s2	2ph - 2ph or vapor
% Cc < 0:  non-wetting_____
%   Ccap>-1		n2	liq. film - meniscus - vapor
%   Ccap<=-1, kkc<=1	n1	liq. film - liq. flow - vapor
%   Ccap<=-1, kkc>1	n3	liq. film - 2ph - vapor
% Cc > 0:  wetting_____
% Cc > :
%   kkc>1		w3	2ph - vapor
%   kkc<=1, kkc>=curv1	w2	liq. flow - vapor
%   kkc<=curv1		w1	liq. film - liq. flow - vapor
% Ccap<= Ccc:  capillary condensation_____
%   kkc<kfkc		k1	liq. film - liq. flow - meniscus
%   kkc<=1 or Ccap<=curv2 k2	liq. flow - meniscus
%   kkc>1, Ccap>curv2	k3	2ph - liq. flow - meniscus
%
% see Fig. 1;

% first the logic;
% Schneider's solutions
if costheta==0
  if kkc<=1
    id = 's1';	% liquid film - liq.flow - vapor or 2ph
  else
    id = 's2';	% 2ph - 2ph or vapor
  end
% My solutions (TL, 2007). The possible, two additional solutions in a small
% range kkc < 1 are not given.
else % Ccap~=0
  if n>1
    error([upper(mfilename) ': The downstream state is in the 2ph-region.\n'...
      'Not supported for a contact angle ~= 90°.']);
  return
  end
  % non-wetting cases
  if Cc < 0
    % eq. (31)
    if Cc>=-1-p1pk/pcap;  id = 'n2'; % eq. (32): liq.film - meniscus - vapor
    else
      if kkc<=1;  id = 'n1';	% eq. (15): liq.film - liq.flow - vapor
      else  id = 'n3';		% eq. (33): liq.film - 2ph - vapor
      end
    end
  % wetting cases
  else  % Cc > 0
    % evaporation within the membrane, no cap. condensation
    if Cc > 1
      if kkc>1;  id = 'w3';	% eq. (23) 2ph - vapor
      else
	% eq. (17)
	if kkc >= n*p12/(p1pk+n*p12+pcap)
	  id = 'w2';		% eq. (18) liq.flow - vapor
	else id = 'w1';		% eq. (15) liq.film - liq.flow - vapor
	end
      end
    % Cc <= 1, capillary condensation
    else
      % eq. (26)
      if kkc < n/(1 + pcap*(pk-n*p1)/(p1pk*(pk-n*p12)))
        id = 'k1';		% eq. (24) liq.film - liq.flow - meniscus
      else
		    % eq. (29)
	if kkc<=1 || Cc<=( kkc+p1pk/pcap ) / ( kkc+p1pk*(kkc-n)/(pcap*(1-n)) );
	  id = 'k2';		% eq. (27) liq.flow - meniscus
	else id = 'k3';		% eq. (30) 2ph - liq.flow - meniscus
	end
      end
    end
  end  % Ccap<0
end  % costheta==0

%% Calculate the solutions for the different types of flow.
% mlin is the mass flux according to the equations reported in TL (2007), m is
% the mass flux expressed below. mlin and m should be the same.

% For 2ph-flows the vapor content is calculated.
% if id == 's2' 'n3' 'w3' 'k3'
if strcmp('3',id(2)) | strcmp(id,'s2')
  % eq. (19) is an implicit eq. for a.
  ares = @(a) 1-fxdot(a) - fnu2ph(a)*fk2ph(a)/(kap*rk*dpk);
  a = fzero(ares,[0 1],optimset('fzero'));
  xdot = fxdot(a); x = fx(a); nu2ph = fnu2ph(a); % k2ph = fk2ph(a);
  % For now, nu2ph is evaluated at T1, p1; This should probably, according to
  % effective vapor viscosity nuvap, improved by using a mean pressure.
  % also for nu2ph apparent vapor viscosity and molecular flow correction should
  % be calculated.
end

switch id
  case 's2'
%   disp('ws: 2ph - vapor');

    % calculate m
    if n>=1
      % downstream state lies in 2ph-region
      n=1;
    end
    nn = 1-n+n*nuvap/nu2ph;
    m = nn*p12*kap/(L*nuvap); mlin = m;
    % de
    ded = n*nuvap/(nu2ph*nn);
    pe = p1 - m*nu2ph*ded*L/kap;

    % write this solution;
    fl.sol.len = -2;
    fl.sol.a3 = a;
    fl.sol.q3 = (1-xdot)*rk*m;
    fl.sol.de = ded*L;
    fl.sol.df = 0;
    fl.lin.deL = ded;
    fl.lin.dfL = 0;
    fl.lin.a3 = a;
    fl.lin.x3 = x;
    [fl.flow(1:2).z] = deal([0 fl.sol.de],[fl.sol.de L]);
    [fl.flow(1:2).T] = deal([T1 T2],[T2 T2]);
    [fl.flow(1:2).p] = deal([p1 pe],[pe p2]);
    [fl.flow(1:2).q] = deal([fl.sol.q3 fl.sol.q3],[0 0]);
    [fl.flow(1:2).a] = deal([a a],[1 1]);
    [fl.flow(1:2).x] = deal([fl.lin.x3 fl.lin.x3],[1 1]);
    [fl.flow(1:2).color] = deal('y','m');

  case 's1'
%   disp('ws: liq.film - liq.flow - vapor or 2ph');

    if n>=1
      disp('    downstream state is a 2ph-mixture');
      n=1;
    end
    nn = 1-n+n*nuvap/nuliq;
    m = nn*p12*kap/(L*nuvap);  mlin = m;
    % de
    ded = n*nuvap/(nuliq*nn);
    dfd = -n*nuvap*kliq*(kapc/kap-1)/(nn*nuliq*kml);
    pe = p1 - m*nuliq*ded*L/kap;

    % write this, unique, solution; return
    %  fl.info.kapf = kapc; however, it stays empty
    fl.sol.len = -3;
    fl.sol.a3 = 0;
    fl.sol.q3 = rk*m;
    fl.sol.de = ded*L;
    fl.sol.df = dfd*L;
    fl.lin.deL = ded;
    fl.lin.dfL = dfd;
    fl.lin.a3 = 0;
    fl.lin.x3 = 0;
    % still a calculation
    T4 = T1+fl.sol.q3*fl.sol.df/kliq;
    [fl.flow(1:3).z]...
      = deal([fl.sol.df 0],[0 fl.sol.de],[fl.sol.de L]);
    [fl.flow(1:3).T] = deal([T1 T4],[T4 T2],[T2 T2]);
    [fl.flow(1:3).p] = deal([p1 p1],[p1 pe],[pe p2]);
    [fl.flow(1:3).q]...
      = deal([fl.sol.q3 fl.sol.q3],[fl.sol.q3 fl.sol.q3],[0 0]);
    [fl.flow(1:3).a] = deal([0 0],[0 0],[1 1]);
    [fl.flow(1:3).x] = deal(fl.flow(1:3).a);
    [fl.flow(1:3).color] = deal('c','c','m');

  case 'w3'
%   disp('2ph - vapor');
    % eq. (23)
    mlin = (kap/L) * ( (pk-p2-n*p12)/nuvap + (p1pk+n*p12)/nu2ph );

    % The second solution for kkc < 1 is neglected here.
    p2ph = p12*n+p1pk;
    pv = p12-p2ph;
    mde = p2ph*kap/nu2ph;
    m = (pv*kap/nuvap+mde)/L;

    % the thermal boundary layer in front of the membrane
    TkT1 = p1pk/dpk;
    % 0.36788=1/e
    Tup = T1+TkT1*[0.1 0.36788 1];
    z1 = - kgas/(m*cpgas);
    zup = [2.3*z1 z1 0];

    % write solution; return
    fl.sol.len = -3;
    fl.sol.p0 = p1;
    fl.sol.a3 = a;
%    fl.sol.q3 = rk*m*(1-xdot);
    fl.sol.de = mde/m;
    fl.sol.df = 0;
    fl.lin.deL = fl.sol.de/L;
    fl.lin.dfL = 0;
    fl.lin.a3 = a;
%    fl.lin.x3 = x;
    [fl.flow(1:3).z] = deal(zup,[0 fl.sol.de],[fl.sol.de L]);
    [fl.flow(1:3).T] = deal(Tup,[Tup(3) T2],[T2 T2]);
    [fl.flow(1:3).p] = deal([p1 p1 p1],[p1 p1-p2ph],[p2+pv p2]);
    [fl.flow(1:3).q] = deal([0 0 0],[fl.sol.q3 fl.sol.q3],[0 0]);
    [fl.flow(1:3).a] = deal([1 1 1],[a a],[1 1]);
    %[fl.flow(1:3).x] = deal([1 1 1],fl.lin.x3*[1 1],[1 1]);
    [fl.flow(1:3).color] = deal('m','y','m');

  case 'w2'
%   disp('liq.flow - vapor');
    % eq. (18)
    mlin=(kap/L) * ( (pk-p2-n*p12)/nuvap + ...
      (pcap+p1pk+n*p12*(1+pcap/p1pk))/((1+kkc*pcap/p1pk)*nuliq) );

    pv = p12*(1-n)-p1pk;
    mde = (pcap+p12*n*(1+pcap/p1pk)+p1pk)...
      / (nuliq/kap+pcap*nuliq/(p1pk*kapc));
    T31 = mde*rk/kml - T12;
    pl = mde*nuliq/kap;
    m = (pv*kap/nuvap+mde)/L;

    % 0.36788=1/e
    Tup = T1+T31*[0.36788 1];
    z1 = - kgas/(m*cpgas);
    zup = [z1 0];

    % write solution; return
    fl.sol.len = -3;
    fl.sol.a3 = 0;
    fl.sol.q3 = rk*m;
    fl.sol.de = mde/m;
    fl.sol.df = 0;
    fl.lin.deL = fl.sol.de/L;
    fl.lin.dfL = 0;
    fl.lin.a3 = 0;
    fl.lin.x3 = 0;
    [fl.flow(1:3).z] = deal(zup,[0 fl.sol.de],[fl.sol.de L]);
    [fl.flow(1:3).T] = deal(Tup,[Tup(2) T2],[T2 T2]);
    [fl.flow(1:3).p] = deal(p1*[1 1],[p2+pv-pcap+pl p2+pv-pcap],[p2+pv p2]);
    [fl.flow(1:3).q] = deal([0 0],[fl.sol.q3 fl.sol.q3],[0 0]);
    [fl.flow(1:3).a] = deal([1 1],[0 0],[1 1]);
    [fl.flow(1:3).x] = deal(fl.flow(1:3).a);
    [fl.flow(1:3).color] = deal('m','c','m');

  case {'w1','n1'}
%   disp('liq.film - liq.flow - vapor');
    % eq. (15)
    mlin=(kap/L) * ( (pk-p2-n*p12)/nuvap + (pcap+p1pk+n*p12)/nuliq );

    pv = p12*(1-n)-p1pk;
    pl = p12+pcap-pv;
    mde = pl*kap/nuliq;
    m = (pv*kap/nuvap+mde)/L;
    T42 = mde*rk/kml;
    df = kliq*(T42-T12)/(m*rk);
    T4 = T42+T2;

    % write solution; return
    fl.sol.len = -3;
    fl.sol.a3 = 0;
    fl.sol.q3 = rk*m;
    fl.sol.de = mde/m;
    fl.sol.df = df;
    fl.lin.deL = fl.sol.de/L;
    fl.lin.dfL = fl.sol.df/L;
    fl.lin.a3 = 0;
    fl.lin.x3 = 0;
    [fl.flow(1:3).z] = ...
      deal([fl.sol.df 0],[0 fl.sol.de],[fl.sol.de L]);
    [fl.flow(1:3).T] = deal([T1 T4],[T4 T2],[T2 T2]);
    [fl.flow(1:3).p] = deal([p1 p1],[p1 p1-pl],[p2+pv p2]);
    [fl.flow(1:3).q] = ...
      deal([fl.sol.q3 fl.sol.q3],[fl.sol.q3 fl.sol.q3],[0 0]);
    [fl.flow(1:3).a] = deal([0 0],[0 0],[1 1]);
    [fl.flow(1:3).x] = deal(fl.flow(1:3).a);
    [fl.flow(1:3).color] = deal('c','c','m');

  case 'n2'
%   disp('liq.film - vapor');
    % eq. (32)
    % here, nuvap is (slightly) wrong
    mlin = (p12*kap/(nuvap*L)) * ( 1 - (n*p1/pk)/(1+(p1pk/pcap)*(1-n*p12/pk)) );

    pv = p12 - n*p12/(1+p1pk/pcap);
    % small error in m:
    m = pv*kap/(nuvap*L);
    df = kliq*T12/(m*rk);

    % write the solution
    fl.sol.len = -2;
    fl.sol.a3 = 0;
    fl.sol.q3 = rk*m;
    fl.sol.de = 0;
    fl.sol.df = -df;
    fl.lin.deL = 0;
    fl.lin.dfL = fl.sol.df/L;
    fl.lin.a3 = 0;
    fl.lin.x3 = 0;
    [fl.flow(1:2).z] = deal([fl.sol.df fl.sol.de],[fl.sol.de L]);
    [fl.flow(1:2).T] = deal([T1 T2],[T2 T2]);
    [fl.flow(1:2).p] = deal([p1 p1],[p2+pv p2]);
    [fl.flow(1:2).q] = deal([fl.sol.q3 fl.sol.q3],[0 0]);
    [fl.flow(1:2).a] = deal([0 0],[1 1]);
    [fl.flow(1:2).x] = deal(fl.flow.a);
    [fl.flow(1:2).color] = deal('c','m');

  case 'n3'
%   disp('liq.film - 2ph - vapor');
    % eq. (33)
    mlin = (kap/L)* ( (n*p12+p1pk)*(1/nu2ph-1/nuvap) + p12/nuvap + pcap/nu2ph );

    T14 = -(p1pk+pcap)/dpk;
    p2ph = dpk*(T12-T14);
    pv = p12+pcap-p2ph;
    mde = p2ph*kap/nu2ph;
    m = (pv*kap/nuvap+mde)/L;
    df = kliq*T14/(m*rk);
    T4 = T1-T14;

    % write solution; return
    fl.sol.len = -3;
    fl.sol.a3 = a;
    fl.sol.q3 = rk*m;
    fl.sol.de = mde/m;
    fl.sol.df = -df;
    fl.lin.deL = fl.sol.de/L;
    fl.lin.dfL = fl.sol.df/L;
    fl.lin.a3 = a;
    fl.lin.x3 = x;
    [fl.flow(1:3).z] = deal([fl.sol.df 0],[0 fl.sol.de],[fl.sol.de L]);
    [fl.flow(1:3).T] = deal([T1 T4],[T4 T2],[T2 T2]);
    [fl.flow(1:3).p] = deal([p1 p1],[p2+pv+p2ph p2+pv],[p2+pv p2]);
    [fl.flow(1:3).q] = deal(fl.sol.q3*[1 1],(1-xdot)*rk*[1 1],[0 0]);
    [fl.flow(1:3).a] = deal([0 0],[a a],[1 1]);
    [fl.flow(1:3).x] = deal([0 0],fl.lin.x3*[1 1],[1 1]);
    [fl.flow(1:3).color] = deal('c','y','m');

  case 'k3'
%   disp('2ph - liq.flow');
    % eq. (30)
    mlin = ( p12*kap/(nuliq*L*(kkc-1)) ) * ( (1-nuliq/nu2ph)*(n+p1pk/p12) ...
     - (1-kkc*nuliq/nu2ph)*(1+pcap*(pk-n*p1)/(p1pk*(pk-n*p12)) - pcap/p12) );

    kck = 1/kkc;
    % the solution
    p2ph = (p12*(1-kck*n+pcap*(1-n)/p1pk)-pcap-kck*p1pk)/(1-kck);
    pl = p12*(1+pcap*(1-n)/p1pk)-p2ph-pcap;
    mdf = p2ph*kap/nu2ph;
    % This also has a small error, not much larger than a rounding error!
    % Expressions for m and mlin checked with mathematica, they are not equal.
    m = (pl*kap/nuliq+mdf)/L;
    T35 = p2ph/dpk;

    % the thermal boundary layer in front of the membrane
    TkT1 = p1pk/dpk;
    % 0.36788=1/e
    Tup = T1+TkT1*[0.1 0.36788 1];
    z1 = - kgas/(m*cpgas);
    zup = [2.3*z1 z1 0];
    T5 = TkT1-T35+T1;

    % write solution; return
    fl.sol.len = -3;
    fl.sol.p0 = p1;
    fl.sol.a3 = a;
    fl.sol.q3 = rk*m*(1-xdot);
    fl.sol.de = L;
    fl.sol.df = mdf/m;
    fl.lin.deL = 1;
    fl.lin.dfL = fl.sol.df/L;
    fl.lin.a3 = a;
    fl.lin.x3 = x;
    [fl.flow(1:3).z] = deal(zup,[0 fl.sol.df],[fl.sol.df L]);
    [fl.flow(1:3).T] = deal(Tup,[Tup(3) T5],[T5 T2]);
    [fl.flow(1:3).p] = deal([p1 p1 p1],[p1 p1-p2ph],...
		[p1-p2ph-pcap p1-p2ph-pcap-p2ph]);
    [fl.flow(1:3).q] = deal([0 0 0],[fl.sol.q3 fl.sol.q3],rk*m*[1 1]);
    [fl.flow(1:3).a] = deal([1 1 1],[a a],[0 0]);
    [fl.flow(1:3).x] = deal([1 1 1],fl.lin.x3*[1 1],[0 0]);
    [fl.flow(1:3).color] = deal('m','y','c');

  case 'k2'
%   disp('liq.flow');
    % eq. (27)
    mlin = (p12*kap/(nuliq*L)) * ( p1pk/pcap+n+(pk-n*p1)/(pk-n*p12) ) ...
      / (kkc+p1pk/pcap);

    m = mlin;
    % This seems erroneous. (By not more than a few percent, however.)
    %m = p12*(1+pcap/p1pk) / (L*nuliq*(pcap/(kapc*p1pk)+1/kap));
    T41 = m*L*rk/kml-T12;
    pl = m*L*nuliq/kap;
    p14 = pcap*dps*T41/p1pk;

    % the thermal boundary layer in front of the membrane
    % 0.36788=1/e
    Tup = T1+T41*[0.1 0.36788 1];
    z1 = - kgas/(m*cpgas);
    zup = [2.3*z1 z1 0];

    % write solution; return
    fl.sol.len = -2;
    fl.sol.p0 = p1;
    fl.sol.a3 = 0;
    fl.sol.q3 = rk*m;
    fl.sol.de = L;
    fl.sol.df = 0;
    fl.lin.deL = 1;
    fl.lin.dfL = 0;
    fl.lin.a3 = 0;
    fl.lin.x3 = 0;
    [fl.flow(1:2).z] = deal(zup,[0 L]);
    [fl.flow(1:2).T] = deal(Tup,[Tup(3) T2]);
    [fl.flow(1:2).p] = deal([p1 p1 p1],[p1-p14 p1-p14-pl]);
    [fl.flow(1:2).q] = deal([0 0 0],rk*m*[1 1]);
    [fl.flow(1:2).a] = deal([1 1 1],[0 0]);
    [fl.flow(1:2).x] = deal([1 1 1],[0 0]);
    [fl.flow(1:2).color] = deal('m','c');

  case 'k1'
%   disp('liq.film - liq.flow');
    % eq. (24)
    mlin = (kap*p12)/(nuliq*L) * ( 1+pcap*(pk-n*p1)/(p1pk*(pk-n*p12)) );
    m = mlin;
    disp([upper(mfilename) ': rudimentary flow struct written, m = mlin.'])
    fl.sol.len = -2; fl.sol.p0 = p1;
    [fl.flow(1:2).color] = deal('c','c');
  otherwise
    error([upper(mfilename) ': The program should never come here.']);
end

fl.info.m = m;
fl.sol.pe = p2;
fl.lin.m = mlin;
fl.info.substance = s; fl.info.membrane = mem; fl.info.fmodel = f;
%------------------------------------------------------------------ end mlinpsat
