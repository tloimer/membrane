function w = stoffwerte(property,T,p)
% STOFFWERTE(property,T,p) returns the value of a thermophysical
% property at given T or T and p.

% check input
if (nargin==2 & (property=='v' | property=='??')
  msg=sprintf(['STOFFWERTE(%s,%g): not enough arguments, '...
    'pressure is missing.'],property,T);
  error(msg);
end

% universale Konstanten
R=8314.4;
% Molekulargewicht nur für v nötig - siehe v.

switch(property)
  case 'rho' % T
    switch(stoff)
      case 'ethanol'
        t0= 1.16239e3;
        t1=-2.25788;
        t2=5.30621e-3;
        t3=-6.63070e-6;
      otherwise notfound;
    end
    w = t0 + t1*T + t2*T^2 + t3*T^3;
  case 'v' % T,p
    switch(stoff)
      case 'ethanol'
        M=46.07;
        b0=9.6838e3;
        b1=-1.3575e7;
        b2=6.3248e9;
        b3=-1.0114e12;
      otherwise notfound;
    end
    B = bo + b1/T + b2/T^2 + b3/T^3;
    a = R*T*1000/p;
    w = (a/2 + sqrt(a^2/4+a*B))/(1000*M);
  case {'ps','dpsdT'}
    switch(stoff)
      case 'ethanol'
        Aa= 6.923365 + 3; % convert to Pa, not kPa
        Ab=1410.46;
        Ac=-64.636;
      otherwise notfound;
    end
    ps = 10^(Aa-Ab/(Ac+T));
    switch(property)
      case 'ps'
        w = ps;
      case 'dpsdT'
        w = ps*Ab*2.302585092994/(Ac+T)^2;
    end
  case 'r' % T,p
    w = T*(v(T,ps(T))-1/rho(T))/stoffwerte('dpsdT',T);
  case {'mul','nul'}
    switch(stoff)
      case 'ethanol'
      otherwise notfound;
  case {'mug','nug'}
  case 'kl'
  case 'kg'
  case 'cpl'
  case 'cpg'
end
 
function notfound
  if (nargin==3)
    msg=sprintf(['STOFFWERTE(%s,%g,%g): Property %s\n'...
      'for substance %s not found.'],property,T,p,property,stoff);
  else
    msg=sprintf(['STOFFWERTE(%s,%g): Property %s\n'...
      'for substance %s not found.'],property,T,,property,stoff);
  end
  error('msg');
