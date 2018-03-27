function [m pK ps] = misotherm(T1,p1,p2,theta,s,mem,ps,pK)
%MISOTHERM  Isothermal mass flux through a homogeneous membrane.
%  MISOTHERM(T1,P1,P2,THETA,S,MEM) returns the mass flux [kg/m2s] of
%  SUBSTANCE S through the homogeneous MEMBRANE MEM with a contact angle
%  THETA according to an isothermal description.
%
%  MISOTHERM(T1,P1,P2,THETA,S,MEM,PS,PK) uses the saturation pressure PS =
%  psat(T1) and the vapor pressure at a curved meniscus PK = pkelv(T1).
%
%  [M PK PS] = MISOTHERM(T1,P1,P2,THETA,S,MEM) returns PK and PS, in
%  addition to the mass flux. If PS and PK are supplied, they are copied to
%  the output unchanged.

% for p1 > pK				    =======================
%					  1 ) 3    liq.   7 ( 8    2
% p1-p2 = p1-p3 + p3-p7 + p7-p8 + p8-p2     =======================
%         `-,-´   `-,-´   `-,-´   `-,-´
% (1)      *       pliq   -2sig/r  pvap    *|	p1-p3 : 2sig/r = ps-p1 : ps-pK
%      p8 = pK => pvap = p8-pK		    |	p1-p3 = (2sig/r)(ps-p1)/(ps-pK)
% (2)  pliq = m nuliq de / kappa
% (3)  pvap = m nuapp (L-de) / kappa
%
% (1,2): pliq = p1-p2 - (2sig/r)*(ps-p1)/(ps-pK) + 2sig/r - (pK-p2)
%  =>  pliq = p1-pK + (2sig/r)(p1-pK)/(ps-pK)
% (2,3): pliq/nuliq + pvap/nuapp = m L / kappa
%	 => m = (kappa/L) ( pliq/nuliq + pvap/nuapp )
%
% for p2 > pK				  ========================
% p1-p2 = p1-p3 + p3-p7 + p7-p2		1 ) 3     liq.         7 ( 8
%         `-,-´   `-,-´   `-,-´		  ========================
% (4)	    *	   pliq     *	   * |  p1-p3 see above
% (5)  pliq = m nuliq L / kappa	   * |  p7-p2 = -(2sig/r)(ps-p2)/(ps-pK)
%
% (4,5): pliq = p1-p2 - (2sig/r)(ps-p1)/(ps-pK) + (2sig/r)(ps-p2)/(ps-pK)
%             = (p1-p2) * ( 1 + (2sig/r)/(ps-pK) )

if nargin == 6
  ps = s.ps(T1);
end

% Some input sanitizing.
if ps < p1
  error([upper(mfilename)...
	': The upstream state is a liquid. This is not implemented.']);
elseif p2 < 0.
  error([upper(mfilename)...
	': The downstream pressure is negative. That is not possible.']);
elseif T1 < 0.
  error([upper(mfilename)...
	': The upstream temperature is below absolute zero. Impossible.']);
end

costheta = cos(theta);
pcap = mem.fcurv(costheta)*s.sigma(T1);
if nargin ~=8
  pK = ps*exp(-pcap/(s.R*s.rho(T1)*T1));
end

p12 = p1 - p2;
betakn_nu = mem.beta * 3*sqrt(pi/(8*s.R))/mem.dia;
kap_L = mem.kappa/mem.L;  corrKn = betakn_nu/sqrt(T1);
if p1 < pK % pure vapor flow
  % m = (p1-p2)*kappa/(nuapp*L)
  m = kap_L*p12*( 1/s.nug(T1,(p1+p2)/2) + corrKn );
else % p1 > pK
  if p2 < pK % & p1 > pK, liq. and vapor flow
    pv_nu = (pK-p2)*( 1/s.nug(T1,(pK+p2)/2) + corrKn );
    pl_nu = (p1 - pK + pcap*(p1-pK)/(ps-pK)) / s.nul(T1);
    m = kap_L*(pv_nu+pl_nu);
  else % p1 > pK & p2 > pK
    m = kap_L*p12*(1+pcap/(ps-pK))/s.nul(T1);
  end
end
