function m = mlinearnew(p1,p2,T1,theta,s,mem,f)
%MLINEAR    Mass flux from linear theory [kg/m2s].
%  MLINEAR(P1,P2,T1,THETA,S,M,F) calculates the mass flux according to linear
%  theory. For P1 = PSAT(T1), the mass flux is calculated as given in JMS07
%  [Loimer, J. Membr. Sci. 301, pp. 107-117, 2007]. For wetted systems and PK -
%  N*(P1-P2) < P1 < PSAT, the mass flux is obtained from an linear interpolation
%  between MGAS and MLINEAR(P1=PSAT), see JMS11 [Loimer, J. Membr, Sci. 383, pp.
%  104-115, 2011]. P1 [Pa] is the upstream pressure, P2 [Pa] the downstream
%  pressure, T1 the upstream temperature, and THETA is the contact angle. at the
%  upstream temperature. The downstream pressure is P2 [Pa], TROOM [°C] is
%  demanded for diagnostic purposes. The substance S, membrane M and two-phase
%  flow model F are provided via structs, see SUBSTANCE, MEMBRANE and FMODEL.
%  The output is a flowstruct FL. The mass flux is reported in FL.LIN.M (or
%  FL.INFO.M) [kg/m2s].
%
%  See also FLOWSTRUCT, SUBSTANCE, MEMBRANE, FMODEL.
%
%  Calls MLINPSAT.

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
