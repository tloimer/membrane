function fl = mkelv(deltap,T1,L)
%MKELV      Mass flux from linear theory [kg/m2s].
%  MKELV(T,DELTAP,L) calculates the mass flux for initial state T and
%  PS(T) for the applied pressure drop DELTAP and membrane thickness L.
%  Returns the solution in a FLOWSTRUCT. Multiple solutions are
%  returned in a FLOWSTRUCT array.
%
%  Calls K, KL, NUG, NUL, PS, R.

% some variables
p1 = ps(T1);
p2 = p1-deltap;
T2 = T1-jt(T1,p1)*deltap;
% modified viscosity
nueff = (nug(T1,p1) + nug(T2,p1-deltap))/2;
% this should have been better, but for liq.film it seems worse
%nueff = (nug(T2,ps(T2))+nug(T2,p1-deltap))/2;
% originally: 
%nueff = nu(T1,1);
kapc = kappac(T1);
kap = kappa;
dps = dpsdT(T1)*(T1-T2);
if costheta==0
  Ca=0;
else
  pcap = curv*sig(T1);
  Ca = dps/pcap; % = dpsdT(T1)*jt(T1,p1)*deltap/pcap(T1)
end
nuliq = nul(T1); % = nu(T1,0)
kliq = kl(T1);
kml = k(T1,0);

% n indicates the downstream state
n = jt(T1,p1)*dpsdT(T1);
% for Ca~=0, also this should be defined differently
if n>=1
  % downstream state lies in 2ph-region
  n=1;
end

% the flowstruct template
% lin.Te is the temperature at the evaporation front
fl = struct('info',struct('kap',kap,'m',[],'Ca',Ca,'kapc',kapc,...
    'kapf',[],'L',L,'T0',T1,'p0',p1,'dp',deltap,'ph',fmodel),...
  'sol',struct('len',[],'a3',[],'q3',[],'T0',T1,'Te',T2,'p0',p1,...
    'pe',p2,'de',[],'df',[]),...
  'flow',struct('z',{},'T',{},'p',{},'q',{},'a',{},'x',{},'color',{}),...
  'lin',struct('m',[],'Te',T2,'T4',[],'deL',[],'dfL',[],'a3',[],'x3',[]));

% first the logic; also some variables which are only needed for some
% cases, are calculated.
if Ca==0
  % schneider's solutions
  if kap>kapc
    % schneider: 2ph - vapor
    id = 'a';
  else
    % schneider: liq.film - liq.flow - vapor
    id = 'b';
  end
else % Ca~=0
  if n>=1
  % downstream state lies in 2ph-region
  error(['The downstream state is in the 2ph-region.\n'...
      'Not supported for a contact angle ~= 90°.'])
  return
  end
  % some variables
  pk = kelv(T1)*p1;
  dpk = dpkdT(T1)*(T1-T2);
  rk = rkelv(T1);
  kapck = nuliq*kml/(dpkdT(T1)*rk);
  if Ca>0
    % critical permeability for wetting liquid
    kapcl = nuliq*kml/(rk*(dpkdT(T1)+(1-kelv(T1))*p1/(T1-T2)));
    %wetting cases
    if kap>kapck
      % id = 'n' as nonsense:  max.meniscus - liq.flow - vapor
      % pressure jump (double meniscus) - 2ph - vapor
      id = 'c';
    else % kap<=kapcl
      kapf = kapc + kap*(1-sqrt(1+4*Ca^2*kapc/kap))/(2*Ca^2);
      fl.info.kapf = kapf;
      %now calculate the modified kappaf, kapfk
      kapold = kappa;
      kapfk = fzero(@kapfkelv,kapf*[.95 1.05],optimset(...
        optimset('fzero'),'TolX',kapf*1e-3),T1,T1-T2,p1,sig(T1),nuliq);
      kappa(kapold);
      if kap>=kapfk
        % liq.flow - vapor
        id = 'd';
      else
        % liq.film - liq.flow - vapor
        id = 'e';
      end
    end
  else % Ca<0, nonwetting
    % calculate the modified critical permeability
    Cak = (dpk+p1-pk)/pcap;
    if Cak>=-1
      % lig.film - vapor
      id = 'f';
      if kap>kapck
        % lig.film - 2ph - vapor, 2nd solution
        id(2) = 'g';
      end
    else %Ca<-1
      if kap>kapck
        % lig.film - 2ph - vapor
        id = 'g';
      else % kap<=kapck
        % liq.film - liq.flow - vapor
        id = 'h';
      end
    end
  end
end

% we may have several solutions
for i=1:length(id)
if i==2
  disp('second solution');
  fl(2)=fl(1);
end
 
switch id(i)
  case 'a'
    % schneider: 2ph - vapor
    disp('2ph - vapor');

    % determine a from  (1-lambda) nu'  k' / (nu k) = kc/k
    kkc = kapc/(kap*nuliq*kml);

    % calculate vapor content a
    a=fzero(@ares,[0 1],optimset('fzero'),T1,kkc);

    % calculate m
    nn = 1-n+n*nueff/nu(T1,a);
    m = nn*deltap*kap/(L*nueff);
    % de
    ded = n*nueff/(nu(T1,a)*nn);
    pe = p1 - m*nu(T1,a)*ded*L/kap;
  
    % write this solution;
    fl(i).sol.len = 2;
    fl(i).sol.a3 = a;
    fl(i).sol.q3 = (1-xdot(T1,a))*r(T1)*m;
    fl(i).sol.de = ded*L;
    fl(i).sol.df = 0;
    fl(i).lin.deL = ded;
    fl(i).lin.dfL = 0;
    fl(i).lin.a3 = a;
    fl(i).lin.x3 = x(T1,a);
    [fl(i).flow(1:2).z] = deal([0 fl(i).sol.de],[fl(i).sol.de L]);
    [fl(i).flow(1:2).T] = deal([T1 T2],[T2 T2]);
    [fl(i).flow(1:2).p] = deal([p1 pe],[pe p2]);
    [fl(i).flow(1:2).q] = deal([fl(i).sol.q3 fl(i).sol.q3],[0 0]);
    [fl(i).flow(1:2).a] = deal([a a],[1 1]);
    [fl(i).flow(1:2).x] = deal([fl(i).lin.x3 fl(i).lin.x3],[1 1]);
    [fl(i).flow(1:2).color] = deal('y','m');

  case 'b'
    % schneider: liq.film - liq.flow - vapor
    disp('liq.film - liq. flow - vapor');

    nn = 1-n+n*nueff/nuliq;
    m = nn*deltap*kap/(L*nueff);
    % de
    ded = n*nueff/(nuliq*nn);
    dfd = -n*nueff*kliq*(kapc/kap-1)/(nn*nuliq*kml);
    pe = p1 - m*nuliq*ded*L/kap;

    % write this, unique, solution; return
    %  fl(i).info.kapf = kapc; however, it stays empty
    fl(i).sol.len = 3;
    fl(i).sol.a3 = 0;
    fl(i).sol.q3 = r(T1)*m;
      % still a calculation
      T4 = T1+fl(i).sol.q3*fl(i).sol.df/kliq;
    fl(i).sol.de = ded*L;
    fl(i).sol.df = dfd*L;
    fl(i).lin.deL = ded;
    fl(i).lin.dfL = dfd;
    fl(i).lin.a3 = 0;
    fl(i).lin.x3 = 0;
    [fl(i).flow(1:3).z]...
      = deal([fl(i).sol.df 0],[0 fl(i).sol.de],[fl(i).sol.de L]);
    [fl(i).flow(1:3).T] = deal([T1 T4],[T4 T2],[T2 T2]);
    [fl(i).flow(1:3).p] = deal([p1 p1],[p1 pe],[pe p2]);
    [fl(i).flow(1:3).q]...
      = deal([fl(i).sol.q3 fl(i).sol.q3],[fl(i).sol.q3 fl(i).sol.q3],[0 0]);
    [fl(i).flow(1:3).a] = deal([0 0],[0 0],[1 1]);
    [fl(i).flow(1:3).x] = deal(fl(i).flow(1:3).a);
    [fl(i).flow(1:3).color] = deal('c','c','m');

  case 'c'
    % (pressure jump), 2ph - vapor

    % here, kkc is reused! see case 'b'.
    kkc = kapck/(kap*nuliq*kml);
    a=fzero(@ares,[0 1],optimset('fzero'),T1,kkc);

    mde = k(T1,a)*(T1-T2)/(rk*(1-xdot(T1,a)));
    p2ph = mde*nu(T1,a)/kap;
    pv = pk-p2-p2ph;
    m = (pv*kap/nueff+mde)/L;

    % write solution; return
    fl(i).sol.len = 2;
    fl(i).sol.p0 = pk;
    fl(i).sol.a3 = a;
    fl(i).sol.q3 = rk*m*(1-xdot(T1,a));
    fl(i).sol.de = mde/m;
    fl(i).sol.df = 0;
    fl(i).lin.deL = fl(i).sol.de/L;
    fl(i).lin.dfL = 0;
    fl(i).lin.a3 = a;
    fl(i).lin.x3 = x(T1,a);
    [fl(i).flow(1:2).z] = deal([0 fl(i).sol.de],[fl(i).sol.de L]);
    [fl(i).flow(1:2).T] = deal([T1 T2],[T2 T2]);
    [fl(i).flow(1:2).p] = deal([pk pk-p2ph],[p2+pv p2]);
    [fl(i).flow(1:2).q] = deal([fl(i).sol.q3 fl(i).sol.q3],[0 0]);
    [fl(i).flow(1:2).a] = deal([a a],[1 1]);
    [fl(i).flow(1:2).x] = deal(fl(i).lin.x3*[1 1],[1 1]);
    [fl(i).flow(1:2).color] = deal('y','m');

%  case 'n' % as nonsense
%    % max.meniscus - liq.flow - vapor
%    disp('max.meniscus - liq. flow - vapor');

%    pl = dpk+p1-pk;
%    mde = pl*kap/nuliq;
%    m = (kap*(deltap-pl)/nueff+mde)/L;
%
%    % write solution; return
%    fl(i).sol.len = 2;
%    fl(i).sol.a3 = 0;
%    fl(i).sol.q3 = rk*m;
%    fl(i).sol.de = mde/m;
%    fl(i).sol.df = 0;
%    fl(i).lin.deL = fl(i).sol.de/L;
%    fl(i).lin.dfL = 0;
%    fl(i).lin.a3 = 0;
%    fl(i).lin.x3 = 0;
%    [fl(i).flow(1:2).z] = deal([0 fl(i).sol.de],[fl(i).sol.de L]);
%    [fl(i).flow(1:2).T] = deal([T1 T2],[T2 T2]);
%    [fl(i).flow(1:2).p] = deal([p1-pcap p1-pcap-pl],[p2+deltap-pl p2]);
%    [fl(i).flow(1:2).q] = deal([fl(i).sol.q3 fl(i).sol.q3],[0 0]);
%    [fl(i).flow(1:2).a] = deal([0 0],[1 1]);
%    [fl(i).flow(1:2).x] = deal(fl(i).flow(1:2).a);
%    [fl(i).flow(1:2).color] = deal('c','m');

  case 'd'
    % liq.flow - vapor
    disp('liq. flow - vapor');
    y = [kap 0 -nuliq 0; 0 kap nueff -nueff*L;...
      0 0 rk 0; 0 1 0 0] \ [0; 0; kml*(T1-T2); deltap+pk-p1-dpk];
    [pl pv mde m] = deal(y(1),y(2),y(3),y(4));

    % write solution; return
    fl(i).sol.len = 2;
    fl(i).sol.a3 = 0;
    fl(i).sol.q3 = rk*m;
    fl(i).sol.de = mde/m;
    fl(i).sol.df = 0;
    fl(i).lin.deL = fl(i).sol.de/L;
    fl(i).lin.dfL = 0;
    fl(i).lin.a3 = 0;
    fl(i).lin.x3 = 0;
    [fl(i).flow(1:2).z] = deal([0 fl(i).sol.de],[fl(i).sol.de L]);
    [fl(i).flow(1:2).T] = deal([T1 T2],[T2 T2]);
    [fl(i).flow(1:2).p] = deal([p2+pv-pcap+pl p2+pv-pcap],[p2+pv p2]);
    [fl(i).flow(1:2).q] = deal([fl(i).sol.q3 fl(i).sol.q3],[0 0]);
    [fl(i).flow(1:2).a] = deal([0 0],[1 1]);
    [fl(i).flow(1:2).x] = deal(fl(i).flow(1:2).a);
    [fl(i).flow(1:2).color] = deal('c','m');

  case {'e','h'}
    % liq.film - liq.flow - vapor
    disp('liq.film - liq. flow - vapor');
    y = [1 1 0 0 0 0; kap 0 -nuliq 0 0 0;...
      0 kap nueff 0 -nueff*L 0; 0 0 0 rk 0 kliq;...
      0 0 rk 0 0 -kml; 0 1 0 0 0 0] \ ...
      [deltap+pcap;0;0;kliq*T1;-kml*T2;deltap+pk-p1-dpk];
    [pl pv mde mdf m T4] = deal(y(1),y(2),y(3),y(4),y(5),y(6));
    % write solution; return
    fl(i).sol.len = 3;
    fl(i).sol.a3 = 0;
    fl(i).sol.q3 = rk*m;
    fl(i).sol.de = mde/m;
    fl(i).sol.df = -mdf/m;
    fl(i).lin.deL = fl(i).sol.de/L;
    fl(i).lin.dfL = fl(i).sol.df/L;
    fl(i).lin.a3 = 0;
    fl(i).lin.x3 = 0;
    [fl(i).flow(1:3).z] = ...
      deal([fl(i).sol.df 0],[0 fl(i).sol.de],[fl(i).sol.de L]);
    [fl(i).flow(1:3).T] = deal([T1 T4],[T4 T2],[T2 T2]);
    [fl(i).flow(1:3).p] = deal([p1 p1],[p1 p1-pl],[p2+pv p2]);
    [fl(i).flow(1:3).q] = ...
      deal([fl(i).sol.q3 fl(i).sol.q3],[fl(i).sol.q3 fl(i).sol.q3],[0 0]);
    [fl(i).flow(1:3).a] = deal([0 0],[0 0],[1 1]);
    [fl(i).flow(1:3).x] = deal(fl(i).flow(1:3).a);
    [fl(i).flow(1:3).color] = deal('c','c','m');

  case 'f'
    % lig.film - vapor
    disp('liq. film - vapor');
    y = [kap 0 -nueff*L; 0 rk 0; 1 0 0] \...
      [0; kliq*(T1-T2); deltap+pk-p1-dpk];
    [pv mdf m] = deal(y(1),y(2),y(3));

    %disp(sprintf('pmen/pcap = %g',(deltap-pv)/(-curv*sig(T1))));
    % write the solution
    fl(i).sol.len = 2;
    fl(i).sol.a3 = 0;
    fl(i).sol.q3 = r(T1)*m;
    fl(i).sol.de = 0;
    fl(i).sol.df = -mdf/m;
    fl(i).lin.deL = 0;
    fl(i).lin.dfL = fl(i).sol.df/L;
    fl(i).lin.a3 = 0;
    fl(i).lin.x3 = 0;
    [fl(i).flow(1:2).z] = deal([fl(i).sol.df fl(i).sol.de],[fl(i).sol.de L]);
    [fl(i).flow(1:2).T] = deal([T1 T2],[T2 T2]);
    [fl(i).flow(1:2).p] = deal([p1 p1],[p2+pv p2]);
    [fl(i).flow(1:2).q] = deal([fl(i).sol.q3 fl(i).sol.q3],[0 0]);
    [fl(i).flow(1:2).a] = deal([0 0],[1 1]);
    [fl(i).flow(1:2).x] = deal(fl(i).flow.a);
    [fl(i).flow(1:2).color] = deal('c','m');

  case 'g'
    % lig.film - 2ph - vapor
    disp('liq. film - 2ph - vapor');

    % here, kkc is reused! see case 'b'.
    kkc = kapck/(kap*nuliq*kml);
    a=fzero(@ares,[0 1],optimset('fzero'),T1,kkc);

    T4 = Tk(p1);

    y = [1 1 0 0 0; kap 0 -nu(T1,a) 0 0; 0 kap nueff 0 -nueff*L; ...
      0 0 0 rk 0; 0 0 rk*(1-xdot(T1,a)) 0 0] \ ...
      [deltap; 0; 0; kliq*(T1-T4); k(T1,a)*(T4-T2)];
    [p2ph pv mde mdf m] = deal(y(1),y(2),y(3),y(4),y(5));

    % write solution; return
    fl(i).sol.len = 3;
    fl(i).sol.a3 = a;
    fl(i).sol.q3 = rk*m;
    fl(i).sol.de = mde/m;
    fl(i).sol.df = -mdf/m;
    fl(i).lin.deL = fl(i).sol.de/L;
    fl(i).lin.dfL = fl(i).sol.df/L;
    fl(i).lin.a3 = a;
    fl(i).lin.x3 = x(T1,a);
    [fl(i).flow(1:3).z] = deal([fl(i).sol.df 0],[0 fl(i).sol.de],[fl(i).sol.de L]);
    [fl(i).flow(1:3).T] = deal([T1 T4],[T4 T2],[T2 T2]);
    [fl(i).flow(1:3).p] = deal([p1 p1],[p2+pv+p2ph p2+pv],[p2+pv p2]);
    [fl(i).flow(1:3).q] = deal(fl(i).sol.q3*[1 1],(1-xdot(T1,a))*rk*[1 1],[0 0]);
    [fl(i).flow(1:3).a] = deal([0 0],[a a],[1 1]);
    [fl(i).flow(1:3).x] = deal([0 0],fl(i).lin.x3*[1 1],[1 1]);
    [fl(i).flow(1:3).color] = deal('c','y','m');

  otherwise
    %error('The program should never come here.');
end

fl(i).info.m = m;
fl(i).sol.pe = p2;
fl(i).lin.m = m;

end %for i=1:length(id)

%-------------------------------------------------------

function ar = ares(a,T1,kkc)
ar = (1-xdot(T1,a)) - nu(T1,a).*k(T1,a).*kkc;
% Schneider:
% xr = (1-x) - nu(T1,x).*k(T1,x).*kkc;

function dp = kapfkelv(kap,T1,dT,p1,sigma,nuliq)
kappa(kap);
dp = dpkdT(T1)*dT+ps(T1)*(1-kelv(T1)) + curv*sigma - ...
  nuliq*k(T1,0)*dT/(kap*rkelv(T1));
