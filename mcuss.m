function fl = mcuss(p1,p12,T1,L)
%MCUSS      Mass flux from isothermal linear theory [kg/m2s].
%  MCUSS(P1,P12,T,L) calculates the mass flux according to an isothermal
%  description (Sidhu & Cussler, J. Mem. Sci., 2001) for the upstream
%  state P1 and T, the applied pressure drop (Delta-) P12 and the
%  membrane thickness L.  Returns the solution in a FLOWSTRUCT.
%
%  Calls K, KL, NUG, NUL, PS, R.

% we set dpkdT=dpsdt and rkelv=r;
% also kapc=kappac, instead of definition with dpkdT and rkelv.

ps1 = ps(T1);
pk1 = kelv(T1)*ps1;
p2 = p1-p12;
kap = kappa;
% try also
%nuvap = (nug(T1,pk1)+nug(T1,p2))/2;
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
  'sol',struct('len',[],'a3',[],'q3',[],'T0',T1,'Te',T1,'p0',p1,...
    'pe',p2,'de',[],'df',[]),...
  'flow',struct('z',{},'T',{},'p',{},'q',{},'a',{},'x',{},'color',{}),...
  'lin',struct('m',[],'Te',T1,'T4',[],'deL',[],'dfL',[],'a3',[],'x3',[]));

if p1<pk1  % vapor flow only

  disp('vapor flow');
  nuvap = (nug(T1,p1)+nug(T1,p2))/2;
  m = p12*kap/(L*nuvap);
  
  % write the result
  fl.sol.len = 1;
  fl.sol.a3 = 1;
  fl.sol.q3 = 0;
  fl.lin.a3 = 1;
  fl.lin.x3 = 1;
  [fl.flow(1).z] = [0 L];
  [fl.flow(1).T] = [T1 T1];
  [fl.flow(1).p] = [p1 p2];
  [fl.flow(1).q] = [0 0];
  [fl.flow(1).a] = [1 1];
  [fl.flow(1).x] = [1 1];
  [fl.flow(1).color] = 'm';

else  % condensation

  disp('liq. flow - vapor');
  pcap = curv*sig(T1);
  nuvap = (nug(T1,pk1)+nug(T1,p2))/2;

  pv = pk1-p2;
  pl = p12-pv+pcap*(1-(ps1-p1)/(ps1-pk1));
  mde = pl*kap/nul(T1);
  m = (pv*kap/nuvap+mde)/L;

  % write solution; return
  fl.sol.len = 2;
  fl.sol.a3 = 0;
  fl.sol.q3 = r(T1)*m;
  fl.sol.de = mde/m;
  fl.lin.deL = fl.sol.de/L;
  fl.lin.a3 = 0;
  fl.lin.x3 = 0;
  [fl.flow(1:2).z] = deal([0 fl.sol.de],[fl.sol.de L]);
%  [fl.flow(1:2).T] = deal([T1 T1],[T1 T1]);
  [fl.flow(1:2).p] = deal([p2+pv+pl-pcap p2+pv-pcap],[p2+pv p2]);
  [fl.flow(1:2).q] = deal([fl.sol.q3 fl.sol.q3],[0 0]);
  [fl.flow(1:2).a] = deal([0 0],[1 1]);
  [fl.flow(1:2).x] = deal(fl.flow(1:2).a);
  [fl.flow(1:2).color] = deal('c','m');

end

fl.sol.df = 0;
fl.lin.dfL = 0;
fl.info.m = m;
fl.lin.m = m;
