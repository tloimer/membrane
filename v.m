function v = v(T,p)
%V(T,p)     Specific volume of the vapor [m3/kg].
%  Calculated by virial equation truncated after the first term.
%
%  Calls B, MOLM.

B=b(T); % [cm3/mol]
[R M]=molm;

v = (R*T./p+B/1000)/M;
