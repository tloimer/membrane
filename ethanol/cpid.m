function cp=cpid(T)
%CPID(T)    Constant pressure heat capacity in ideal gas state [J/kgK].
%
%  Called from CPG.
%
%  Ethanol.
%  Valid for 200K < T < 1500K.
%  From Perry (1997), Table 2-198.

%Compared with data taken from NIST: Thermodynamics Research Center (1997);
%Gurvich et al. (1989).
% t = [50;100;150;200;300;400];
% cm= [37.12;41.7;46.94;52.02;65.49;81.22]*1e3/M;
% cp = interp1q(t,cm,T')';

[R M]=molm;

c1=0.492e5;
c2=1.4577e5;
c3=1.6628e3;
c4=0.939e5;
c5=744.7;

c3t=c3/T;
c5t=c5/T;

cp = ( c1 + c2*(c3t/sinh(c3t))^2 + c4*(c5t/cosh(c5t))^2 )/M;
