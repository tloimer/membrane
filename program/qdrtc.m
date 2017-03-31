function  [Z1, Z2] = qdrtc(A, B, C)
%  [z1, z2] = qdrtc(a, b, c)  solves the quadratic equation
%  a*z^2 - 2*b*z + c == 0  for roots  z1  and  z2  computed
%  accurately out to the last sig. bit or two.  Real arrays
%  of coefficients  a, b, c  yield arrays of roots  z1, z2 .
%  To ease root-tracing,  real roots are ordered:  Z1 >= Z2 .
%  NO PRECAUTIONS AGAINST PREMATURE OVER/UNDERFLOW,  NOR NANS
%  NOR INFINITIES AMONG THE COEFFICIENTS,  NOR  a == 0 ,
%  HAVE BEEN TAKEN IN THIS VERSION OF THE PROGRAM.
%                                        (C) 2004  W. Kahan
%  From https://people.eecs.berkeley.edu/~wkahan/Qdrtcs.pdf,
%  see kahan04.pdf.
sA = size(A) ;  sB = size(B) ;  sC = size(C) ;
sZ = max([sA; sB; sC]) ;  Z1 = ones(sZ(1), sZ(2)) ;
if any( [sZ-sA, sZ-sB, sZ-sC] ) %...  mix scalars and arrays
    if     (sum(sA)==2),  A = A(Z1) ;
    elseif (sum(sB)==2),  B = B(Z1) ;
    elseif (sum(sC)==2),  C = C(Z1) ;
    else  error('Sizes of  qdrtc(A,B,C)''s  arguments mismatch.')
  end,  end
A = A(:) ;  B = B(:) ;  C = C(:) ;  Z1 = Z1(:) ;
if any(imag([A; B; C]))
    error('qdrtc(A, B, C)  accepts only real  A, B, C.'),  end
Z2 = Z1 ;  %...  Allocate initial memory for roots.
D = dscrmt(A, B, C) ;  %...  Discriminant:  see file  dscrmt.m
nD = (D <= 0) ;
if any(nD)  %... Complex conjugate or coincident real roots
    Z1(nD) = B(nD)./A(nD) + sqrt(D(nD))./A(nD) ;
    Z2(nD) = conj(Z1(nD)) ;  end
nD = ~nD ;  if any(nD)  %...  Distinct real roots
    S = B(nD) ;
    S = sqrt(D(nD)).*( sign(S) + (S==0) ) + S ;
    Z1(nD) = S./A(nD) ;  Z2(nD) = C(nD)./S ;  end
nD = (Z1 < Z2) ;  if any(nD)  %...  Sort real roots
    S = Z1(nD) ;  Z1(nD) = Z2(nD) ;  Z2(nD) = S ;  end
Z1 = reshape(Z1, sZ(1), sZ(2)) ;
Z2 = reshape(Z2, sZ(1), sZ(2)) ;
return  %...  End  qdrtc.m
end  %...  End  qdrtc.m

function  D = dscrmt(A, B, C)
%  dscrmt(A, B, C) = B.*B - A.*C  computed extra-precisely
%  if necessary to ensure accuracy out to the last sig. bit
%  or two.  Real columns  A, B, C  must have the same size.
%  This program is intended to work with versions  3.5 - 6.5
%  of  MATLAB  on  PCs,  Power Macs,  and old  680x0 Macs.
%                                     (C) 2004  W. Kahan
%  Determine  Matlab's  Arithmetic Style  AS :
[AS, rl] = arithstyle;

%  Is the obvious way to compute  D  adequately accurate?
if  (AS == 2)  %...  Sum scalar products to  64  sig. bits
    pie = 1024 ;  n = length(A) ;  D = zeros(n,1) ;
    for  j = 1:n
        D(j) = [B(j), A(j)]*[B(j); -C(j)] ;
      end  %...  of loop on  j
  else  %...  All arithmetic rounds to  53  sig. bits
    pie = 3 ;  D = B.*B - A.*C ;  end
E = B.*B + A.*C ;
k = ( pie*abs(D) >= E ) ;
if all(k),  return, end  %... If the obvious way was good enough.
%  Recompute those values of  D  too inaccurate the obvious way.
k = ~k ;
a = A(k) ;  b = B(k) ;  c = C(k) ;
p = b.*b ;  q = a.*c ;  n = length(a) ;
dp = p ;  dq = q ;  %... allocate memory.
if  (AS > 2)  %...  Use the hardware's Fused Multiply-Add
    if rl,  for  j = 1:n
       dp(j) = [b(j), p(j)]*[b(j); -1] ;
   dq(j) = [a(j), q(j)]*[c(j); -1] ; end
      d = (p-q) + (dp-dq) ;  end
  else  %...  Break operands into half-width fragments
    [ah, at] = break2(a) ;  [bh, bt] = break2(b) ;
    [ch, ct] = break2(c) ;
    if  (AS < 2)  %...  All arithmetic rounds to  53  sig. bits
        dp = ((bh.*bh - p) + 2*bh.*bt) + bt.*bt ;
        dq = ((ah.*ch - q) + (ah.*ct + at.*ch)) + at.*ct ;
        d = (p-q) + (dp-dq) ;
      else  %...  Arithmetic may round to  64  and then  53 s.b.
        if rl  %...   Sum scalar products right-to-left
            for j = 1:n
                dp(j) = [bt(j), 2*bh(j), p(j), bh(j)]* ...
                        [bt(j);  bt(j);   -1,  bh(j)] ;
                dq(j) = [at(j), at(j), ah(j), q(j), ah(j)]* ...
                        [ct(j); ch(j); ct(j);  -1;  ch(j)] ;
              end %...  of loop on  j
            d = [dq, dp, q, p]*[-1; 1; -1; 1] ;
          else  %...  Sum scalar products left-to-right
            for j = 1:n
                dp(j) = [bh(j), p(j), 2*bh(j), bt(j)]* ...
                        [bh(j);  -1;   bt(j);  bt(j)] ;
                dq(j) = [ah(j), q(j), ah(j), at(j), at(j)]* ...
                        [ch(j);  -1;  ct(j); ch(j); ct(j)] ;
              end %...  of loop on  j
            d = [p, q, dp, dq]*[1; -1; 1; -1] ;
          end  %...  of extra-precisely summed scalar products
      end  %...  of arithmetic with half-width fragments
  end  %...  Now  d  is fairly accurate.
D(k) = d ;
return  %... End  dscrmt.m
end  %... End  dscrmt.m

function  [Xh, Xt] = break2(X)
%  [Xh, Xt] = break2(X)  produces  Xh = X rounded to 26 sig. bits
%  and  Xt = X - Xh  exactly in 26 sig. bits,  so products like
%  Xh.*Xh,  Xh.*Xt,  Xt.*Xt  can all be computed exactly.  But if
%  arithmetic double-rounds to 64 and then 53 sig. bits,  some of
%  Xt  or  Xh  may require 27 sig. bits for all I know.   W. K.
bigX = X*134217729 ;  %... = X*(2^27 + 1)
Y = X - bigX ;  Xh = Y + bigX ;  %... DON'T OPTIMIZE  Y  AWAY!
Xt = X - Xh ;  return  %... End break2
end   %... End break2


function [arithstyle, rightleft] = arithstyle
persistent AS; persistent rl;
% Once known, bypass the test.
if isempty(AS) || isempty(rl)
%  Determine  Matlab's  Arithmetic Style  AS :
y = 0.5^31 ;  z = -y*y ;
x = [1+y, 1]*[1-y; -1] ;  y = [1, 1+y]*[-1; 1-y] ;
x0 = (x==0) ;  xz = (x==z) ;  y0 = (y==0) ;  yz = (y==z) ;
AS = x0*y0+ 2*xz*yz + 3*x0*yz + 4*xz*y0;
%  Determine whether  MATLAB  adds scalar products right-to-left :
rl = ( [eps, 9, -9]*[eps; 9; 9] > 0 );
if  (AS == 0) | ( (AS == 3) & rl ) | ( (AS == 4) & (~rl) )
  ArithmeticStyle = AS ,  RightToLeft = rl ,
  disp(' Something strange has happened!  Please inform  W. Kahan')
  disp(' about your computer and your version of  Matlab  because')
  error(' dscrmt  did not recognize  Matlab''s  arithmetic style.')
end
end
arithstyle = AS;
rightleft = rl;
end % end arithstyle
