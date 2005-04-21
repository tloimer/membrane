function fl = mlin(deltap,T0,L)
%MLIN       Mass flux from linear theory [kg/m2s].
%  MLIN(DELTAP,T,L) calculates the mass flux for initial state T and
%  ps(T) for the applied pressure drop DELTAP and membrane thickness L.
%  Returns the solution in a FLOWSTRUCT. Multiple solutions are
%  returned in a FLOWSTRUCT array.
%
%  Calls K, KL, NUG, NUL, PS, R.

% some variables
p0 = ps(T0);
Te = T0-jt(T0,p0)*deltap;
% modified viscosity
nueff = (nug(T0,p0) + nug(Te,p0-deltap))/2;
% for liq.film this seems worse
%nueff = (nug(Te,ps(Te))+nug(Te,p0-deltap))/2;
% originally: 
%nueff = nu(T0,1);
kapc = kappac(T0);
kapf = [];
kap = kappa;
dps = dpsdT(T0).*(T0-Te);
if costheta==0
  Ca=[];
else
  Ca = dps/(curv*sig(T0)); % = dpsdT(T0)*jt(T0,p0)*deltap/pcap(T0)
end
nuliq = nul(T0); % = nu(T0,0)
kliq = kl(T0);
kml = k(T0,0);
i = 1; % flowstruct counter

% n indicates the downstream state
n = jt(T0,p0).*dpsdT(T0);
if n>=1
  % downstream state lies in 2ph-region
  n=1;
end

% the flowstruct template
% lin.Te is the temperature at the evaporation front
temp = struct('info',struct('kap',kap,'m',[],'Ca',Ca,'kapc',kapc,...
    'kapf',[],'L',L,'T0',T0,'p0',p0,'dp',deltap,'ph',fmodel),...
  'sol',struct('len',[],'a3',[],'q3',[],'T0',T0,'Te',Te,'p0',p0,...
    'pe',p0-deltap,'de',[],'df',[]),...
  'flow',struct('z',{},'T',{},'p',{},'q',{},'a',{},'x',{},'color',{}),...
  'lin',struct('m',[],'Te',Te,'T4',[],'deL',[],'dfL',[],'a3',[],'x3',[]));

% Logic:
% all n:
%  kap>kapc,  wetting or nonwetting                          2ph - vapor
% only n < 1 (downstream state unsaturated vapor)
%  kap<=kapc,  Ca=0                            liq.film - liq.flow - vapor
%  all kap,   Ca>=-1, nonwetting               lig.film - meniscus - vapor
%  kap<=kapc,  kap>=kapf, wetting   meniscus - liq.flow - meniscus - vapor
%  (kap<=kapc), kap<kapf, wetting   liq.film - liq.flow - meniscus - vapor
%  kap<=kapc,  Ca<-1, nonwetting    liq.film - liq.flow - meniscus - vapor
%

if Ca<0 & Ca>=-1 & n<1
  % ----------------- liq. film - meniscus - vapor
  disp('liq. film - meniscus - vapor');
  y = [kap 0 -nueff*L; 0 r(T0) 0; 1 0 0] \...
    [0; kliq*(T0-Te); deltap-dps];
  [pv mdf m] = deal(y(1),y(2),y(3));
  % write the solution
  fl = temp;
  fl.info.m = m;
  fl.sol.pe = fl(1).sol.pe;
  fl.sol.len = 2;
  fl.sol.a3 = 0;
  fl.sol.q3 = r(T0)*m;
  fl.sol.de = 0;
  fl.sol.df = -mdf/m;
  fl.lin.m = m;
  fl.lin.T4 = Te;
  fl.lin.deL = 0;
  fl.lin.dfL = fl.sol.df/L;
  fl.lin.a3 = 0;
  fl.lin.x3 = 0;
  [fl.flow(1:2).z] = deal([fl.sol.df fl.sol.de],[fl.sol.de L]);
  [fl.flow(1:2).T] = deal([T0 Te],[Te Te]);
  [fl.flow(1:2).p] = deal([p0 p0],[fl.sol.pe+pv fl.sol.pe]);
  [fl.flow(1:2).q] = deal([fl.sol.q3 fl.sol.q3],[0 0]);
  [fl.flow(1:2).a] = deal([0 0],[1 1]);
  [fl.flow(1:2).x] = deal(fl.flow.a);
  [fl.flow(1:2).color] = deal('c','m');
  i = 2;
end

if kap>kapc
  disp('2ph - vapor');
  % ------------ 2ph - vapor
  % supercritical permeability
  % see Schneider (1983)

  % determine a from  (1-lambda) nu'  k' / (nu k) = kc/k
  kkc = kapc./(kap*nuliq.*kml);

  % calculate vapor content a
  a=fzero(@ares,[0 1],optimset('fzero'),T0,kkc);

  % calculate m
  nn = 1-n+n.*nueff./nu(T0,a);
  m = nn.*deltap.*kap./(L.*nueff);
  % de
  ded = n.*nueff./(nu(T0,a).*nn);
  pe = p0 - m.*nu(T0,a).*ded*L/kap;
  
  % write this solution;
  fl(i) = temp;
  fl(i).info.m = m;
  fl(i).sol.len = 2;
  fl(i).sol.a3 = a;
  fl(i).sol.q3 = (1-xdot(T0,a))*r(T0)*m;
  fl(i).sol.de = ded*L;
  fl(i).sol.df = 0;
  fl(i).lin.m = m;
  fl(i).lin.T4 = T0;
  fl(i).lin.deL = ded;
  fl(i).lin.dfL = 0;
  fl(i).lin.a3 = a;
  fl(i).lin.x3 = x(T0,a);
  [fl(i).flow(1:2).z] = deal([0 fl(i).sol.de],[fl(i).sol.de L]);
  [fl(i).flow(1:2).T] = deal([T0 Te],[Te Te]);
  [fl(i).flow(1:2).p] = deal([p0 pe],[pe fl(i).sol.pe]);
  [fl(i).flow(1:2).q] = deal([fl(i).sol.q3 fl(i).sol.q3],[0 0]);
  [fl(i).flow(1:2).a] = deal([a a],[1 1]);
  [fl(i).flow(1:2).x] = deal([fl(i).lin.x3 fl(i).lin.x3],[1 1]);
  [fl(i).flow(1:2).color] = deal('y','m');

else % kappa<=kappac
  % subcritical permeability
  if costheta==0
    disp('liq.film - liq. flow - vapor');
    % ------------ liq.film - liq. flow - vapor
    % see Schneider (1983)
    nn = 1-n+n.*nueff./nuliq;
    m = nn.*deltap.*kap./(L.*nueff);
    % de
    ded = n.*nueff./(nuliq.*nn);
    dfd = -n.*nueff.*kliq.*(kapc/kap-1)./(nn.*nuliq.*kml);
    pe = p0 - m.*nuliq.*ded*L/kap;

    % write this, unique, solution; return
    fl(i) = temp;
    fl(i).info.m = m;
    %  fl(i).info.kapf = kapc; however, it stays empty
    fl(i).sol.len = 3;
    fl(i).sol.a3 = 0;
    fl(i).sol.q3 = r(T0)*m;
    fl(i).sol.de = ded*L;
    fl(i).sol.df = dfd*L;
    fl(i).lin.m = m;
    fl(i).lin.T4 = T0+fl(i).sol.q3*fl(i).sol.df/kliq;
    fl(i).lin.deL = ded;
    fl(i).lin.dfL = -dfd;
    fl(i).lin.a3 = 0;
    fl(i).lin.x3 = 0;
    [fl(i).flow(1:3).z]...
      = deal([fl(i).sol.df 0],[0 fl(i).sol.de],[fl(i).sol.de L]);
    [fl(i).flow(1:3).T] = deal([T0 fl(i).lin.T4],[fl(i).lin.T4 Te],[Te Te]);
    [fl(i).flow(1:3).p] = deal([p0 p0],[p0 pe],[pe fl(i).sol.pe]);
    [fl(i).flow(1:3).q]...
      = deal([fl(i).sol.q3 fl(i).sol.q3],[fl(i).sol.q3 fl(i).sol.q3],[0 0]);
    [fl(i).flow(1:3).a] = deal([0 0],[0 0],[1 1]);
    [fl(i).flow(1:3).x] = deal(fl(i).flow(1:3).a);
    [fl(i).flow(1:3).color] = deal('c','c','m');

  elseif n<=1 % costheta~=0; I only provide for unsaturated state 2
    if Ca>0 % costheta>0, wetting
      kapf = kapc + kap*(1-sqrt(1+4*Ca^2*kapc/kap))/(2*Ca^2);
      if kap>=kapf % just for wetting case
        disp('meniscus - liq. flow - vapor');
        % ------------- meniscus - liq. flow - meniscus - vapor
        y = [kap 0 -nuliq 0; 0 kap nueff -nueff*L;...
          0 0 r(T0) 0; 0 1 0 0] \ [0; 0; kml*(T0-Te); deltap-dps];
        [pl pv mde m] = deal(y(1),y(2),y(3),y(4));

        % write solution; return
        fl(i) = temp;
        fl(i).info.m = m;
        fl(i).info.kapf = kapf;
        fl(i).sol.len = 2;
        fl(i).sol.a3 = 0;
        fl(i).sol.q3 = r(T0)*m;
        fl(i).sol.de = mde/m;
        fl(i).sol.df = 0;
        fl(i).lin.m = m;
        fl(i).lin.T4 = Te;
        fl(i).lin.deL = fl(i).sol.de/L;
        fl(i).lin.dfL = fl(i).sol.df/L;
        fl(i).lin.a3 = 0;
        fl(i).lin.x3 = 0;
        [fl(i).flow(1:2).z] = deal([0 fl(i).sol.de],[fl(i).sol.de L]);
        [fl(i).flow(1:2).T] = deal([T0 Te],[Te Te]);
        [fl(i).flow(1:2).p] = deal([fl(i).sol.pe+pv+pl-curv*sig(T0) p0-curv*sig(T0)-pl],...
	  [fl(i).sol.pe+pv fl(i).sol.pe]);
        [fl(i).flow(1:2).q] = deal([fl(i).sol.q3 fl(i).sol.q3],[0 0]);
        [fl(i).flow(1:2).a] = deal([0 0],[1 1]);
        [fl(i).flow(1:2).x] = deal(fl(i).flow(1:2).a);
        [fl(i).flow(1:2).color] = deal('c','m');
	return % dirty programming - now Ca > 0 means (Ca>0 & kap<kapf)
      end
    end
    if Ca<-1 | Ca>0 % eig: Ca<-1 | ( Ca>0 & kap<kapf )
      disp('liq.film - liq. flow - meniscus - vapor');
      % ------------- liq. film - liq. flow - meniscus - vapor
      y = [1 1 0 0 0 0; kap 0 -nuliq 0 0 0;...
        0 kap nueff 0 -nueff*L 0; 0 0 0 r(T0) 0 kliq;...
	0 0 r(T0) 0 0 -kml; 0 1 0 0 0 0] \ ...
	[deltap+curv*sig(T0);0;0;kliq*T0;-kml*Te;deltap-dps];
      [pl pv mde mdf m T4] = deal(y(1),y(2),y(3),y(4),y(5),y(6));
      % write solution; return
      fl(i) = temp;
      fl(i).info.m = m;
      fl(i).info.kapf = kapf;
      fl(i).sol.len = 3;
      fl(i).sol.a3 = 0;
      fl(i).sol.q3 = r(T0)*m;
      fl(i).sol.de = mde/m;
      fl(i).sol.df = -mdf/m;
      fl(i).lin.m = m;
      fl(i).lin.T4 = T4;
      fl(i).lin.deL = fl(i).sol.de/L;
      fl(i).lin.dfL = fl(i).sol.df/L;
      fl(i).lin.a3 = 0;
      fl(i).lin.x3 = 0;
      [fl(i).flow(1:3).z]...
        = deal([fl(i).sol.df 0],[0 fl(i).sol.de],[fl(i).sol.de L]);
      [fl(i).flow(1:3).T] = deal([T0 T4],[T4 Te],[Te Te]);
      [fl(i).flow(1:3).p]...
        = deal([p0 p0],[p0 p0-pl],[fl(i).sol.pe+pv fl(i).sol.pe]);
      [fl(i).flow(1:3).q]...
        = deal([fl(i).sol.q3 fl(i).sol.q3],[fl(i).sol.q3 fl(i).sol.q3],[0 0]);
      [fl(i).flow(1:3).a] = deal([0 0],[0 0],[1 1]);
      [fl(i).flow(1:3).x] = deal(fl(i).flow(1:3).a);
      [fl(i).flow(1:3).color] = deal('c','c','m');
    end
    % case Ca>=-1, kap<kapc unresolved; was treated first
  else % n>1
    error(['The downstream state is in the 2ph-region.\n'...
      'Not supported for a contact angle ~=0.'])
  end
end

function ar = ares(a,T0,kkc)
ar = (1-xdot(T0,a)) - nu(T0,a).*k(T0,a).*kkc;
% Schneider:
% xr = (1-x) - nu(T0,x).*k(T0,x).*kkc;
