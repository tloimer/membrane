function cp=cpg(T,p)
%CPG(T,P)   Constant pressure heat capacity of the vapor [J/kgK].
%  CPG = CPID + P*DCPDP(T).
%
%  Calls CPID, DCPDP.

%cpold = cpid(T)
dcp = dcpdp(T);
cp = cpid(T) + p.*dcp;

%%%
%  Ethanol
%  Data taken from NIST: Thermodynamics Research Center (1997);
%  Gurvich et al. (1989).

% values given in J/(molK) at 1 bar (?).

%[R M]=molm;

%t = [50;100;150;200;300;400];
%cm= [37.12;41.7;46.94;52.02;65.49;81.22]*1e3/M;

%cp = interp1q(t,cm,T')';
