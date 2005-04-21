function fl = munsat(p1,p12,T1,L)
%MUNSAT     Mass flux from isothermal linear theory [kg/m2s].
%  MUNSAT(P1,P12,T,L) calculates the mass flux for the upstream state P1
%  and T, the applied pressure drop (Delta-) P12 and the membrane
%  thickness L.  Returns the solution in a FLOWSTRUCT.
%
%  Calls K, KL, NUG, NUL, PS, R.

% we set dpkdT=dpsdt and rkelv=r;
% also kapc=kappac, instead of definition with dpkdT and rkelv.

T12 = jt(T1,p1)*p12; % = T1-T2 = mu_jt*p12
T2 = T1 - T12;
Ts1 = Ts(p1);
%ps1 = ps(T1);
p1pk = (1-kelv(Ts1))*p1; % = p1 - pk
dpk = dpsdT(Ts1);
Tk1 = Ts1+p1pk/dpk;
p2 = p1-p12;
kap = kappa;
% for comparison, we use the simplest viscosity
% try also
% nuvap = (nug(T2,p2)+nug(T1,p2))/2;
% modified viscosity
% nuvap = (nug(T1,p1) + nug(T2,p1-p12))/2;
% this should have been better, but for liq.film it seems worse
%nuvap = (nug(T2,ps(T2))+nug(T2,p1-p12))/2;
% originally: 
%nuvap = nu(T1,1);

% the flowstruct template
% lin.Te is the temperature at the evaporation front
fl = struct('info',struct('kap',kap,'m',[],'Ca',[],'kapc',kappac(T1),...
    'kapf',[],'L',L,'T0',T1,'p0',p1,'dp',p12,'ph',fmodel),...
  'sol',struct('len',[],'a3',[],'q3',[],'T0',T1,'Te',T2,'p0',p1,...
    'pe',p2,'de',[],'df',[]),...
  'flow',struct('z',{},'T',{},'p',{},'q',{},'a',{},'x',{},'color',{}),...
  'lin',struct('m',[],'Te',T2,'T4',[],'deL',[],'dfL',[],'a3',[],'x3',[]));


if T1-T12>Tk1  % T2>Tk1, vapor flow only

  disp('vapor flow');
  nuvap = (nug(T2,p1)+nug(T2,p2))/2;
  m = p12*kap/(L*nuvap);

    % the thermal boundary layer in front of the membrane
    % 0.36788=1/e; 0.1=e^(-2.3)
    Tup = T1-T12*[0.1 0.36788 0.7408 1];
    z1 = - kg(T1)/(m*cpg(T1,p1));
    zup = z1*[2.3 1 0.3 0];

  % write the result
% has to be corrected!
  fl.sol.len = 2;
  fl.sol.a3 = 1;
  fl.sol.q3 = 0;
  fl.lin.a3 = 1;
  fl.lin.x3 = 1;
  [fl.flow(1:2).z] = deal(zup,[0 L]);
  [fl.flow(1:2).T] = deal(Tup,[T2 T2]);
  [fl.flow(1:2).p] = deal(p1*[1 1 1 1],[p1 p2]);
  [fl.flow(1:2).q] = deal([0 0 0 0],[0 0]);
  [fl.flow(1:2).a] = deal([1 1 1 1],[1 1]);
  [fl.flow(1:2).x] = deal(fl.flow(1:2).a);
  [fl.flow(1:2).color] = deal('m','m');

else  % condensation
  if kap>kappac(Ts1) error('kappa must be smaller than kappac.');
  end

  disp('liq. flow - vapor');
  pcap = curv*sig(T1);
  TsT2 = T12-T1+Ts1;
  kapc = kappac(Ts1);
  nuliq = nul(Ts1);

  pv = p12-dpk*TsT2-p1pk;
  mde=(pcap+p1pk+dpk*TsT2*(pcap/p1pk+1))/(nuliq*(1/kap+pcap/(kapc*p1pk)));
  T31 = mde*nuliq/(kapc*dpk) - T12;
  pl = mde*nuliq/kap;
  % now the vapor viscosity
  nuvap = (nug(T2,p2+pv)+nug(T2,p2))/2;
  m = (pv*kap/nuvap+mde)/L;

  % the thermal boundary layer in front of the membrane
  % 0.36788=1/e; 0.1=e^(-2.3)
  Tup = T1+T31*[0.1 0.36788 0.7408 1];
  z1 = - kg(T1)/(m*cpg(T1,p1));
  zup = z1*[2.3 1 0.3 0];

  % write solution; return
  fl.sol.len = 3;
  fl.sol.a3 = 0;
  fl.sol.q3 = r(T1)*m;
  fl.sol.de = mde/m;
  fl.lin.deL = fl.sol.de/L;
  fl.lin.a3 = 0;
  fl.lin.x3 = 0;
  [fl.flow(1:3).z] = deal(zup,[0 fl.sol.de],[fl.sol.de L]);
  [fl.flow(1:3).T] = deal(Tup,[T1+T31 T2],[T2 T2]);
  [fl.flow(1:3).p] = deal([p1 p1 p1 p1],[p2+pv+pl-pcap p2+pv-pcap],[p2+pv p2]);
  [fl.flow(1:3).q] = deal([0 0 0 0],[fl.sol.q3 fl.sol.q3],[0 0]);
  [fl.flow(1:3).a] = deal([1 1 1 1],[0 0],[1 1]);
  [fl.flow(1:3).x] = deal(fl.flow(1:3).a);
  [fl.flow(1:3).color] = deal('m','c','m');

end

fl.sol.df = 0;
fl.lin.dfL = 0;
fl.info.m = m;
fl.lin.m = m;
