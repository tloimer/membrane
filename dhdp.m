function dhdp = dhdp(T)
%DHDP(T)    Derivative dh/dp for T constant [m3/kg].
%  DHDP(T) is based on the virial equation truncated after the first
%  term.
%
%  Calls TDBDT -> VIRIAL, MOLM.
%  See also V.

[B TdBdT] = TdbdT(T);
[R M] = molm;

dhdp = (B-TdBdT)/(1000*M);
