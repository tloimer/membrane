function ds = dsigdT(T)
%DSIGDT(T)  Derivative of surface tension [N/mK].
%
%  A polynom of order 2 was fit to the data. The derivative of this
%  polnomial fit is taken.
%
%  Ethanol.
%  Valid for 340K < T < 415K.
%  From Landolt-Börnstein: Group IV, vol. 16 (1997).
%
%  See SIG.

%sig = 3.8028e-2 - 2.2643e-5*T - 1.0275e-7*T.^2;
ds =  - 2.2643e-5 - 2*1.0275e-7*T;
