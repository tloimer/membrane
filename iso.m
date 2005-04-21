function fl = mkelv(p12)
%MKELV      Mass flux from linear theory [kg/m2s].
%  MKELV(P12,T,L) calculates the mass flux for the upstream state T and
%  PS(T), the applied pressure drop (Delta-) P12 and the membrane
%  thickness L.  Returns the solution in a FLOWSTRUCT.
%
%  Calls K, KL, NUG, NUL, PS, R.

% we set dpkdT=dpsdt and rkelv=r;
% also kapc=kappac, instead of definition with dpkdT and rkelv.

% some variables
L=25e-6;
T1=293;
p1 = 3e5; %ps(T1);
p2 = p1-p12;
T12 = 2.8e-5*p12; % = T1-T2 = mu_jt*p12
T2 = T1-T12;
kapc = 1.7e-17; %kappac(T1);
% should be:
% kapc = nu(T,0).*k(T,0)./(dpkdT(T).*rkelv(T));
kap = 1.6e-17;
kkc = kap/kapc;
%dps = dpsdT(T1)*T12;

% n indicates the downstream state
n = 0.26;%jt(T1,p1)*dpsdT(T1);

p1pk = 0.126e5; %(1-kelv(T1))*p1; % = p1 - pk

% now the vapor viscosity
nuvap = 4.3e-7; %(nug(T2,p1-n*p12-p1pk)+nug(T2,p2))/2;
% p1-n*p12 = p2+p12*(n-1)
%nuvap = (nug(T1,p1) + nug(T2,p1-p12))/2;
% this should have been better, but for liq.film it seems worse
%nuvap = (nug(T2,ps(T2))+nug(T2,p1-p12))/2;
% originally: 
%nuvap = nu(T1,1);

if costheta~=0
  pcap = 1e6; %curv*sig(T1);
  Ccap = (n*p12 + p1pk)/pcap;
end


% the flowstruct template
% lin.Te is the temperature at the evaporation front
fl = struct('info',struct('kap',kap,'m',[],'Ca',Ccap,'kapc',kapc,...
    'kapf',[],'L',L,'T0',T1,'p0',p1,'dp',p12,'ph',fmodel),...
  'sol',struct('len',[],'a3',[],'q3',[],'T0',T1,'Te',T2,'p0',p1,...
    'pe',p2,'de',[],'df',[]),...
  'flow',struct('z',{},'T',{},'p',{},'q',{},'a',{},'x',{},'color',{}),...
  'lin',struct('m',[],'Te',T2,'T4',[],'deL',[],'dfL',[],'a3',[],'x3',[]));

% the logic:
% costheta==0:  Schneider's solutions_____
%   kkc<=1		s1	liquid film - liq. flow - vapor or 2ph
%   kkc>1		s2	2ph - 2ph or vapor
% Ccap < 0:  non-wetting_____
%   Ccap>-1		n2	liq. film - meniscus - vapor 
%   Ccap<=-1, kkc<=1	n1	liq. film - liq. flow - vapor
%   Ccap<=-1, kkc>1	n3	liq. film - 2ph - vapor
% Ccap > 0:  wetting_____
% Ccap > Ccc:
%   kkc>1		w3	2ph - vapor
%   kkc<=1, kkc>=curv1	w2	liq. flow - vapor
%   kkc<=curv1		w1	liq. film - liq. flow - vapor
% Ccap<= Ccc:  capillary condensation_____
%   kkc<kfkc		k1	liq. film - liq. flow - meniscus
%   kkc<=1 or Ccap<=curv2 k2	liq. flow - meniscus
%   kkc>1, Ccap>curv2	k3	2ph - liq. flow - meniscus

% first the logic; also some variables which are only needed for some
% cases, are calculated.
if costheta==0
  if kkc<=1
    id = 's1';
  else
    id = 's2';
  end
else % Ccap~=0
  if n>1  error(['The downstream state is in the 2ph-region.\n'...
      'Not supported for a contact angle ~= 90°.'])
  return
  end

  if Ccap < 0  % non-wetting cases______
    if Ccap>-1  id = 'n2';
    else
      if kkc<=1  id = 'n1';
      else  id = 'n3';
      end
    end
  else  % Ccap > 0, wetting______
    Ccc = p1pk/(pcap*(1-n));
    if Ccap > Ccc  % wetting flow______
      if kkc>1  id = 'w3';
      else
        if kkc >= (Ccap-p1pk/pcap)/(1+Ccap)
	  id = 'w2';
	else id = 'w1';
	end
      end
    else  % Ccap<=Ccc capillary condensation______
      if kkc < n/(1+1/Ccc)
        id = 'k1';
      else
	if kkc<=1 | Ccap<=((n/kkc+1)*p1pk/pcap+1)/(1+1/Ccc)
	  id = 'k2';
	else id = 'k3';
	end
      end
    end
  end  % Ccap<0
end  % costheta==0

% some abbreviations
nuliq = 3.2e-7; %nul(T1); % = nu(T1,0)
kliq = 0.104;
kml = 0.17; %k(T1,0);
 
switch id
  case 's2'
    % schneider: 2ph - vapor
    disp('ws: 2ph - vapor');

    % determine a from  (1-lambda) nu'  k' / (nu k) = kc/k
    rdpk_k = kapc/(kap*nuliq*kml);
    % calculate vapor content a
    a=fzero(@ares,[0 1],optimset('fzero'),T1,rdpk_k);

    % calculate m
    if n>=1
      % downstream state lies in 2ph-region
      n=1;
    end
    nn = 1-n+n*nuvap/nu(T1,a);
    m = nn*p12*kap/(L*nuvap);
    % de
    ded = n*nuvap/(nu(T1,a)*nn);
    pe = p1 - m*nu(T1,a)*ded*L/kap;
  
    % write this solution;
    fl.sol.len = 2;
    fl.sol.a3 = a;
    fl.sol.q3 = (1-xdot(T1,a))*r(T1)*m;
    fl.sol.de = ded*L;
    fl.sol.df = 0;
    fl.lin.deL = ded;
    fl.lin.dfL = 0;
    fl.lin.a3 = a;
    fl.lin.x3 = x(T1,a);
    [fl.flow(1:2).z] = deal([0 fl.sol.de],[fl.sol.de L]);
    [fl.flow(1:2).T] = deal([T1 T2],[T2 T2]);
    [fl.flow(1:2).p] = deal([p1 pe],[pe p2]);
    [fl.flow(1:2).q] = deal([fl.sol.q3 fl.sol.q3],[0 0]);
    [fl.flow(1:2).a] = deal([a a],[1 1]);
    [fl.flow(1:2).x] = deal([fl.lin.x3 fl.lin.x3],[1 1]);
    [fl.flow(1:2).color] = deal('y','m');

  case 's1'
    % schneider: liq.film - liq.flow - vapor
    disp('ws: liq.film - liq. flow - vapor or 2ph');

    if n>=1
      % downstream state lies in 2ph-region
      n=1;
    end
    nn = 1-n+n*nuvap/nuliq;
    m = nn*p12*kap/(L*nuvap);
    % de
    ded = n*nuvap/(nuliq*nn);
    dfd = -n*nuvap*kliq*(kapc/kap-1)/(nn*nuliq*kml);
    pe = p1 - m*nuliq*ded*L/kap;

    % write this, unique, solution; return
    %  fl.info.kapf = kapc; however, it stays empty
    fl.sol.len = 3;
    fl.sol.a3 = 0;
    fl.sol.q3 = r(T1)*m;
      % still a calculation
      T4 = T1+fl.sol.q3*fl.sol.df/kliq;
    fl.sol.de = ded*L;
    fl.sol.df = dfd*L;
    fl.lin.deL = ded;
    fl.lin.dfL = dfd;
    fl.lin.a3 = 0;
    fl.lin.x3 = 0;
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
    disp('2ph - vapor');

    % calculate the vapor content
    rdpk_k = kapc/(kap*nuliq*kml);
    a=fzero(@ares,[0 1],optimset('fzero'),T1,rdpk_k);

    % the solution
    p2ph = p12*n+p1pk;
    mde = p2ph*kap/nu(T1,a);
    pv = p12-p2ph;
    % now the vapor viscosity
    nuvap = (nug(T2,p2+pv)+nug(T2,p2))/2;
    m = (pv*kap/nuvap+mde)/L;

    % the thermal boundary layer in front of the membrane
    TkT1 = p1pk/dpsdT(T1);
    % 0.36788=1/e
    Tup = T1+TkT1*[0.1 0.36788 1];
    z1 = - kg(T1)/(m*cpg(T1,p1));
    zup = [2.3*z1 z1 0];

    % write solution; return
    fl.sol.len = 3;
    fl.sol.p0 = p1;
    fl.sol.a3 = a;
    fl.sol.q3 = r(T1)*m*(1-xdot(T1,a));
    fl.sol.de = mde/m;
    fl.sol.df = 0;
    fl.lin.deL = fl.sol.de/L;
    fl.lin.dfL = 0;
    fl.lin.a3 = a;
    fl.lin.x3 = x(T1,a);
    [fl.flow(1:3).z] = deal(zup,[0 fl.sol.de],[fl.sol.de L]);
    [fl.flow(1:3).T] = deal(Tup,[Tup(3) T2],[T2 T2]);
    [fl.flow(1:3).p] = deal([p1 p1 p1],[p1 p1-p2ph],[p2+pv p2]);
    [fl.flow(1:3).q] = deal([0 0 0],[fl.sol.q3 fl.sol.q3],[0 0]);
    [fl.flow(1:3).a] = deal([1 1 1],[a a],[1 1]);
    [fl.flow(1:3).x] = deal([1 1 1],fl.lin.x3*[1 1],[1 1]);
    [fl.flow(1:3).color] = deal('m','y','m');

  case 'w2'
    disp('liq. flow - vapor');
    %rk = r(T1);
    rk = 340e3;

    pv = p12*(1-n)-p1pk;
    mde = (pcap+p12*n*(1+pcap/p1pk)+p1pk)...
      / (nuliq/kap+pcap*nuliq/(p1pk*kapc));
    T31 = mde*rk/kml - T12;
    pl = mde*nuliq/kap;
    % now the vapor viscosity
    nuvap = (nug(T2,p2+pv)+nug(T2,p2))/2;
    m = (pv*kap/nuvap+mde)/L;

    % 0.36788=1/e
    Tup = T1+T31*[0.36788 1];
    z1 = - kg(T1)/(m*cpg(T1,p1));
    zup = [z1 0];

    % write solution; return
    fl.sol.len = 3;
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
    disp('liq.film - liq. flow - vapor');
    rk = r(T1);

    pv = p12*(1-n)-p1pk;
    pl = p12+pcap-pv;
    mde = pl*kap/nuliq;
    % now the vapor viscosity
    nuvap = (nug(T2,p2+pv)+nug(T2,p2))/2;
    m = (pv*kap/nuvap+mde)/L;
    T42 = mde*rk/kml;
    df = kliq*(T42-T12)/(m*rk);

    T4 = T42+T2;

    % write solution; return
    fl.sol.len = 3;
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
    disp('liq. film - vapor');
    rk = r(T1);

    pv = p12 - n*p12/(1+p1pk/pcap);
    % now the vapor viscosity
    nuvap = (nug(T2,p2+pv)+nug(T2,p2))/2;
    m = pv*kap/(nuvap*L);
    df = kliq*T12/(m*rk);

    % write the solution
    fl.sol.len = 2;
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
    disp('liq. film - 2ph - vapor');
    dpk = dpsdT(T1);
    rk = r(T1);

    rdpk_k = kapc/(kap*nuliq*kml);
    a=fzero(@ares,[0 1],optimset('fzero'),T1,rdpk_k);

    T14 = -(p1pk+pcap)/dpsdT(T1);
    p2ph = dpk*(T12-T14);
    pv = p12+pcap-p2ph;
    mde = p2ph*kap/nu(T1,a);
    % now the vapor viscosity
    nuvap = (nug(T2,p2+pv)+nug(T2,p2))/2;
    m = (pv*kap/nuvap+mde)/L;
    df = kliq*T14/(m*rk);

    T4 = T1-T14;

    % write solution; return
    fl.sol.len = 3;
    fl.sol.a3 = a;
    fl.sol.q3 = rk*m;
    fl.sol.de = mde/m;
    fl.sol.df = -df;
    fl.lin.deL = fl.sol.de/L;
    fl.lin.dfL = fl.sol.df/L;
    fl.lin.a3 = a;
    fl.lin.x3 = x(T1,a);
    [fl.flow(1:3).z] = deal([fl.sol.df 0],[0 fl.sol.de],[fl.sol.de L]);
    [fl.flow(1:3).T] = deal([T1 T4],[T4 T2],[T2 T2]);
    [fl.flow(1:3).p] = deal([p1 p1],[p2+pv+p2ph p2+pv],[p2+pv p2]);
    [fl.flow(1:3).q] = deal(fl.sol.q3*[1 1],(1-xdot(T1,a))*rk*[1 1],[0 0]);
    [fl.flow(1:3).a] = deal([0 0],[a a],[1 1]);
    [fl.flow(1:3).x] = deal([0 0],fl.lin.x3*[1 1],[1 1]);
    [fl.flow(1:3).color] = deal('c','y','m');

  case 'k3'
    disp('2ph - liq. flow');
    dpk = dpsdT(T1);
    rk = r(T1);

    % calculate the vapor content
    rdpk_k = kapc/(kap*nuliq*kml);
    a=fzero(@ares,[0 1],optimset('fzero'),T1,rdpk_k);

    kck = 1/kkc;
    % the solution
    p2ph = (p12*(1-kck*n+pcap*(1-n)/p1pk)-pcap-kck*p1pk)/(1-kck);
    pl = p12*(1+pcap*(1-n)/p1pk)-p2ph-pcap;
    mdf = p2ph*kap/nu(T1,a);
    m = (pl*kap/nuliq+mdf)/L;
    T35 = p2ph/dpk;

    % the thermal boundary layer in front of the membrane
    TkT1 = p1pk/dpk;
    % 0.36788=1/e
    Tup = T1+TkT1*[0.1 0.36788 1];
    z1 = - kg(T1)/(m*cpg(T1,p1));
    zup = [2.3*z1 z1 0];

    T5 = TkT1-T35+T1;

    % write solution; return
    fl.sol.len = 3;
    fl.sol.p0 = p1;
    fl.sol.a3 = a;
    fl.sol.q3 = rk*m*(1-xdot(T1,a));
    fl.sol.de = L;
    fl.sol.df = mdf/m;
    fl.lin.deL = 1;
    fl.lin.dfL = fl.sol.df/L;
    fl.lin.a3 = a;
    fl.lin.x3 = x(T1,a);
    [fl.flow(1:3).z] = deal(zup,[0 fl.sol.df],[fl.sol.df L]);
    [fl.flow(1:3).T] = deal(Tup,[Tup(3) T5],[T5 T2]);
    [fl.flow(1:3).p] = deal([p1 p1 p1],[p1 p1-p2ph],...
		[p1-p2ph-pcap p1-p2ph-pcap-p2ph]);
    [fl.flow(1:3).q] = deal([0 0 0],[fl.sol.q3 fl.sol.q3],rk*m*[1 1]);
    [fl.flow(1:3).a] = deal([1 1 1],[a a],[0 0]);
    [fl.flow(1:3).x] = deal([1 1 1],fl.lin.x3*[1 1],[0 0]);
    [fl.flow(1:3).color] = deal('m','y','c');

  case 'k2'
    disp('liq. flow');
    rk = r(T1);

    m = p12*(1+pcap/p1pk) / (L*nuliq*(pcap/(kapc*p1pk)+1/kap));
    T41 = m*L*rk/kml-T12;
    pl = m*L*nuliq/kap;
    p14 = pcap*dpsdT(T1)*T41/p1pk;

    % the thermal boundary layer in front of the membrane
    % 0.36788=1/e
    Tup = T1+T41*[0.1 0.36788 1];
    z1 = - kg(T1)/(m*cpg(T1,p1));
    zup = [2.3*z1 z1 0];

    % write solution; return
    fl.sol.len = 2;
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

  otherwise
    error('The program should never come here.');
end

fl.info.m = m;
fl.sol.pe = p2;
fl.lin.m = m;

%-------------------------------------------------------

function ar = ares(a,T1,rdpk_k)
ar = (1-xdot(T1,a)) - nu(T1,a).*k(T1,a).*rdpk_k;
% Schneider:
% xr = (1-x) - nu(T1,x).*k(T1,x).*rdpk_k;
