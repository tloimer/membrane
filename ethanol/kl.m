function kl = kl(T)
%KL(T)      Thermal conductivity of the liquid [W/mK].
%
%  Ethanol.
%  Valid for 223.15K < T < 450K.
%  From Perry (1997), Table 2-256 and Table 2-370.
%  Perry also gives an estimation technique, p2-368

%  Data merged at [363.15,410].

%  From Perry (1997), Table 2-256.
%t= [  250;  280;  290;  300;  330;  400;  410;  430;  440;  450];
%k =[0.177;0.171;0.170;0.168;0.159;0.145;0.144;0.140;0.139;0.137];

%  From Perry (1997), Table 2-370.
%t=[223.15;243.15;253.15;303.15;363.15;373.15];
%k=[ 0.188; 0.184; 0.181; 0.171; 0.153;0.151];


t=[223.15;243.15;253.15;303.15;363.15;  410;  430;  440;  450];
k=[ 0.188; 0.184; 0.181; 0.171; 0.153;0.144;0.140;0.139;0.137];
kl = interp1q(t,k,T')';
