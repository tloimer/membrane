function kg = kg(T)
%KG(T)      Thermal conductivity of the vapor [W/mK].
%  =:-0  spooky!  =:-0
%
%  Ethanol.
%  Valid for 250K < T < 400K.
%  Data taken from CRC Handbook.
%  Data is extrapolated (250) and linear (not prop sqrt(T) )
%  interpolation is used.
%  Perry (1997) gives an estimation technique.

t = [250;300;400];
k= [8.7;14.4;25.8]*1e-3;

kg = interp1q(t,k,T')';
