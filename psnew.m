function ps = psnew(T)
%PSNEW(T)      Vapor pressure [Pa].
%
%  Ethanol.
%  Valid for 159 K < T < 513.92 K.
%  From Perry (1997), Table 2-6.
%  P [Pa] = exp( c1 + c2/T + c3*ln(T) + c4*T^c5 ), [T] = K.

ps = exp(74.475-7164.3./T-7.327*log(T)+3.134e-6*T.^2);
