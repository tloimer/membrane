function J = thomas2001(p1,p2,T,M,eta,K0,B0,L,r1)
%THOMAS2001 Molar flux through layered plane or cylindrical membranes.
%  THOMAS2001(P1,P2,T,M,ETA,K0,B0,L) returns the molar flux [kmol/m2]
%  of an ideal gas through a plane ceramic membrane according to the
%  integration of equation (2) in Thomas et al. (2001).
%
%  THOMAS2001(P1,P2,T,M,ETA,K0,B0,L,R1) returns the molar flow per unit
%  length [kmol/m] through a cylindrical membrane. The inner radius
%  is R1, and the first element in K0, B0 or L refers to the first,
%  innermost layer.
%
%  K0, B0, and L must all be scalars, or vectors of the same length.
%
%  Plenty of (derived) results are printed.
%
%  See Thomas et al., Catal. Today 67, 205-216 (2001), or
%  equation (2) in Hussain et al., J. Porous Media 12, 749-757 (2009).

% Notation:
%   J   molar flux, used for plane geometry
%   n   molar flow per unit length, in cylindrical geometry
%
% The molar flux is given by
%
%           1     4       /8 R T |   B0      dp
%    J = - ---  ( - K0   / ----- | + --- p ) --  .
%          R T    3    \/  pi M      eta     dz
%
% Let us introduce two abbreviations, a, b, such that
%
%   J  = - (a + b p) dp/dz.

if nargin > 8, iscyl = true; else iscyl = false; end;
R = 8314.4;	% J/kmolK
RT = R*T;
a = 4.*K0*sqrt(8./(pi*M*RT))/3.;
b = B0./(eta*RT);

% Integration is straightforward in plane geometry, where J = const,
%  
%   J (z2 - z1) = a (p1 - p2) + b ((p1 + p2)/2) (p1 - p2).              (1)
%
% For cylindrical geometry, the molar flux depends on r, but the molar
% flow per unit length, n, is constant.
%
%   J = n / (2 pi r),
%
% integration yields
%
%   (n/(2 pi)) ln(r2/r1) =  a (p1 - p2) + b ((p1 + p2)/2) (p1 - p2).    (2)
%
% For purely molecular flow, for plane geometry where the length of
% the ith layer (a_i, b_i) is given by
%
%   L_i = z_(i+1) - z_i,
%
%      n              n
%   J sum(L_i/a_i) = sum(p_i-p_(i+1)) = p_1 - p_n.
%     i=1            i=1

Jknudsen = (p1 - p2) / sum(L./a);

% For cylindrical geometry
%
%    n    n     r_i+1          n
%   ---- sum(ln(-----)/a_i) = sum(p_i-p_(i+1)) = p_1 - p_n.
%   2 pi i=1     r_i          i=1

if iscyl
    ro = r1 + cumsum(L);
    if size(ro,1) == 1
	ri = [r1 ro(1:end-1)];
    else
	ri = [r1; ro(1:end-1)];
    end
    nknudsen = 2*pi*(p1 - p2) / sum(reallog(ro./ri)./a);
end

% For purely viscous flow, in plane geometry
%
%        n              n
%   2*J sum(L_i/b_i) = sum(p_i^2 - p_(i+1)^2) = p_1^2 - p_n^2.
%       i=1            i=1
%
% while in cylindrical geometry,
%
%   n   n     r_i+1          n
%   -- sum(ln(-----)/b_i) = sum(p_i^2 - p_(i+1)^2) = p_1^2 - p_n^2.
%   pi i=1     r_i          i=1
%

Jviscous = 0.5*(p1*p1 - p2*p2) / sum(L./b);
if iscyl
    nviscous = pi*(p1*p1 - p2*p2) / sum(reallog(ro./ri)./b);
end

% Because addition is commutative, the flow direction does not matter
% for purely viscous or purely molecular (Knudsen) flow.

% v = R*T/(p*M), J = kmol/m2 ( * kg/kmol * m3/kg )
STP = R*273.15/1e5;
fprintf('STP: volume at p = 10^5 Pa, T = 273.15 K.\n');
fprintf('Purely molecular flow,\n');
fprintf('  plane geometry:       J/delta p = %g m^3(STP)/h/bar/m^2\n',...
	Jknudsen*STP*3600/(p1 - p2));
if iscyl
    fprintf(['  cylindrical geometry: J_dp(ri) = %g, '...
	    'J_dp(ro) = %g m^3(STP)/h/bar/m^2\n',...
	    nknudsen*STP*3600/((2*pi*r1)*(p1-p2)),...
	    nknudsen*STP*3600/((2*pi*ro(end))*(p1-p2)));
end
fprintf('\nPurely viscous flow,\n');
fprintf('  plane geometry:       J/delta p = %g m^3(STP)/h/bar/m^2\n',...
	Jviscous*STP*3600/(p1 - p2));
if iscyl
    fprintf(['  cylindrical geometry: J_dp(ri) = %g, '...
	    'J_dp(ro) = %g m^3(STP)/h/bar/m^2\n',...
	    nviscous*STP*3600/((2*pi*r1)*(p1-p2)),...
	    nviscous*STP*3600/((2*pi*ro(end))*(p1-p2)));
end

% The combination of viscous and molecular flow is solved by finding the
% zero of a function, that computes the pressures, given J or n.
% For plane geometry, from eq. (1) above, with pii = p2 and pi = p1,
%   J (z2 - z1) = a (p1 - p2) + b ((p1 + p2)/2) (p1 - p2).          (1)
%
%   J*Li = ai*pi - ai*pii + (bi/2)*pi^2 - (bi/2)*pii^2,
%
%           ai      /ai^2   2                 |
%   pii = - -- +   / ---- + -- (ai*pi - J*Li) |  .
%           bi   \/  bi^2   bi
%
% For cylindrial geometry, from  eq.(2) follows
%   (n/(2 pi)) ln(r2/r1) =  a (p1 - p2) + b ((p1 + p2)/2) (p1 - p2).    (2)
function pii = pres(J)
    pa = p1;
    for i = 1:length(L)
	t = a(i)/b(i);
	pii = -t + sqrt(t*t + 2.*(a(i)*pa - J*L(i))/b(i));
	pa = pii;
    end
end

J = fzero(pres,Jviscous+Jknudsen);

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END THOMAS2001 %%%
