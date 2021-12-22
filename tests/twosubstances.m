function twosubstances(name1, name2)
%TWOSUBSTANCES Compare substance NAME1 and NAME2, with respect to nist.
%  TWOSUBSTANCES(NAME1, NAME2) creates a number of figures which show the
%  difference in the values of properties for two substances.
%
%  See also TESTSUBSTANCE

% Exported form http://webbook.nist.gov/chemistry/fluid to tab delimited data.
% Columns for saturation properties - temperature increments:
% 1 Temperature (K),  2 Pressure (MPa),  3 Density (l, kg/m3),
% 4 Volume (l, m3/kg),  5 Internal Energy (l, kJ/kg),  6 Enthalpy (l, kJ/kg),
% 7 Entropy (l, J/g*K),  8 Cv (l, J/g*K),  9 Cp (l, J/g*K),
% 10 Sound Spd. (l, m/s),  11 Joule-Thomson (l, K/MPa),  12 Viscosity (l, Pa*s),
% 13 Therm. Cond. (l, W/m*K),  14 Surf. Tension (l, N/m), 15 Density (v, kg/m3),
% 16 Volume (v, m3/kg),  17 Internal Energy (v, kJ/kg),  18 Enthalpy (v, kJ/kg),
% 19 Entropy (v, J/g*K),  20 Cv (v, J/g*K),  21 Cp (v, J/g*K),
% 22 Sound Spd. (v, m/s),  23 Joule-Thomson (v, K/MPa),  24 Viscosity (v, Pa*s),
% 25 Therm. Cond. (v, W/m*K)


switch name1
case {'propane', 'propaneperry'}

% Propane, propaneperry, from commit 87546e1.
% Data retrieved on 30 Sep. and 20. Okt. 2014 from the National Institute of
% Standards website, http://webbook.nist.gov/chemistry/fluid.

%Tc = 369.82;
% Saturation properties - temperature increments.
% Data retrieved 2014-09-30
% Simply export data to tab delimited data
% column numbers see above
%
%Temperature (K)	Pressure (MPa)	Density (l, kg/m3)	Volume (l, m3/kg)	Internal Energy (l, kJ/kg)	Enthalpy (l, kJ/kg)	Entropy (l, J/g*K)	Cv (l, J/g*K)	Cp (l, J/g*K)	Sound Spd. (l, m/s)	Joule-Thomson (l, K/MPa)	Viscosity (l, Pa*s)	Therm. Cond. (l, W/m*K)	Surf. Tension (l, N/m)	Density (v, kg/m3)	Volume (v, m3/kg)	Internal Energy (v, kJ/kg)	Enthalpy (v, kJ/kg)	Entropy (v, J/g*K)	Cv (v, J/g*K)	Cp (v, J/g*K)	Sound Spd. (v, m/s)	Joule-Thomson (v, K/MPa)	Viscosity (v, Pa*s)	Therm. Cond. (v, W/m*K); ...
satdata = ...
[100	2.5330e-08	718.20	0.0013924	-168.55	-168.55	-1.0910	1.3238	1.9012	2030.4	-0.62852	0.0037803	0.20323	0.035435	1.3434e-06	7.4440e+05	360.27	379.13	4.3858	0.74794	0.93650	153.65	420.53	2.9792e-06	0.0024171; ...
110	3.4728e-07	708.03	0.0014124	-149.43	-149.43	-0.90878	1.3219	1.9233	1964.8	-0.61839	0.0022603	0.19916	0.033874	1.6744e-05	59724.	367.94	388.68	3.9831	0.78463	0.97319	160.39	308.20	3.2205e-06	0.0029439; ...
120	2.9622e-06	697.87	0.0014329	-130.09	-130.09	-0.74056	1.3229	1.9434	1896.2	-0.60874	0.0015032	0.19451	0.032320	0.00013092	7638.3	375.96	398.59	3.6651	0.81971	1.0083	166.83	233.46	3.4667e-06	0.0034998; ...
130	1.7600e-05	687.73	0.0014540	-110.56	-110.56	-0.58426	1.3253	1.9622	1827.2	-0.59904	0.0010805	0.18940	0.030775	0.00071805	1392.7	384.33	408.84	3.4111	0.85337	1.0420	172.99	181.79	3.7167e-06	0.0040849; ...
140	7.8928e-05	677.59	0.0014758	-90.847	-90.847	-0.43816	1.3290	1.9808	1758.6	-0.58882	0.00082234	0.18394	0.029239	0.0029904	334.41	393.02	419.41	3.2065	0.88602	1.0747	178.91	144.88	3.9700e-06	0.0046989; ...
150	0.00028325	667.43	0.0014983	-70.943	-70.943	-0.30085	1.3342	2.0001	1690.7	-0.57766	0.00065284	0.17821	0.027713	0.010019	99.813	402.02	430.29	3.0407	0.91821	1.1071	184.60	117.77	4.2258e-06	0.0053418; ...
160	0.00084980	657.22	0.0015216	-50.841	-50.839	-0.17112	1.3412	2.0208	1623.4	-0.56521	0.00053462	0.17230	0.026198	0.028195	35.467	411.32	441.46	2.9057	0.95045	1.1398	190.03	97.395	4.4832e-06	0.0060131; ...
170	0.0022035	646.94	0.0015457	-30.521	-30.518	-0.047937	1.3502	2.0436	1556.7	-0.55115	0.00044801	0.16625	0.024694	0.068881	14.518	420.90	452.89	2.7956	0.98328	1.1736	195.20	81.782	4.7415e-06	0.0067123; ...
180	0.0050673	636.56	0.0015709	-9.9619	-9.9539	0.069570	1.3615	2.0690	1490.5	-0.53518	0.00038195	0.16012	0.023203	0.14988	6.6720	430.73	464.54	2.7057	1.0171	1.2090	200.07	69.615	4.9999e-06	0.0074387; ...
190	0.010547	626.06	0.0015973	10.864	10.881	0.18217	1.3751	2.0974	1424.7	-0.51702	0.00032994	0.15397	0.021725	0.29641	3.3737	440.80	476.39	2.6322	1.0524	1.2467	204.61	59.995	5.2578e-06	0.0081917; ...
200	0.020194	615.40	0.0016250	31.988	32.021	0.29051	1.3912	2.1291	1359.3	-0.49635	0.00028793	0.14782	0.020262	0.54151	1.8467	451.08	488.37	2.5723	1.0894	1.2875	208.74	52.290	5.5148e-06	0.0089710; ...
210	0.036037	604.55	0.0016541	53.443	53.503	0.39519	1.4097	2.1646	1294.1	-0.47283	0.00025330	0.14172	0.018814	0.92597	1.0800	461.53	500.45	2.5235	1.1285	1.3318	212.43	46.045	5.7708e-06	0.0097770; ...
220	0.060583	593.47	0.0016850	75.266	75.368	0.49671	1.4307	2.2042	1229.1	-0.44601	0.00022428	0.13572	0.017383	1.4981	0.66752	472.13	512.57	2.4840	1.1697	1.3803	215.59	40.931	6.0262e-06	0.010611; ...
230	0.096787	582.12	0.0017179	97.494	97.660	0.59552	1.4541	2.2481	1164.3	-0.41538	0.00019964	0.12981	0.015970	2.3136	0.43224	482.83	524.67	2.4521	1.2133	1.4337	218.19	36.704	6.2821e-06	0.011477; ...
240	0.14801	570.45	0.0017530	120.17	120.43	0.69203	1.4798	2.2970	1099.5	-0.38021	0.00017847	0.12403	0.014578	3.4357	0.29106	493.60	536.68	2.4264	1.2594	1.4927	220.14	33.184	6.5401e-06	0.012378; ...
250	0.21798	558.40	0.0017908	143.33	143.72	0.78661	1.5077	2.3513	1034.8	-0.33958	0.00016010	0.11840	0.013206	4.9365	0.20257	504.40	548.56	2.4060	1.3080	1.5580	221.39	30.235	6.8027e-06	0.013323; ...
260	0.31070	545.90	0.0018318	167.03	167.60	0.87960	1.5379	2.4118	970.07	-0.29217	0.00014401	0.11293	0.011859	6.8985	0.14496	515.18	560.22	2.3897	1.3592	1.6309	221.88	27.757	7.0729e-06	0.014321; ...
270	0.43046	532.88	0.0018766	191.32	192.13	0.97132	1.5702	2.4796	905.14	-0.23616	0.00012979	0.10763	0.010537	9.4181	0.10618	525.89	571.60	2.3768	1.4129	1.7128	221.54	25.672	7.3552e-06	0.015389; ...
280	0.58173	519.23	0.0019259	216.26	217.38	1.0621	1.6046	2.5563	839.90	-0.16892	0.00011711	0.10252	0.0092437	12.610	0.079299	536.47	582.60	2.3664	1.4694	1.8060	220.30	23.923	7.6551e-06	0.016544; ...
290	0.76921	504.82	0.0019809	241.91	243.43	1.1522	1.6411	2.6446	774.19	-0.086512	0.00010569	0.097597	0.0079822	16.618	0.060175	546.82	593.11	2.3580	1.5287	1.9137	218.06	22.463	7.9800e-06	0.017814; ...
300	0.99780	489.48	0.0020430	268.36	270.40	1.2421	1.6799	2.7482	707.76	0.017115	9.5302e-05	0.092868	0.0067563	21.626	0.046242	556.86	603.00	2.3508	1.5905	2.0406	214.72	21.268	8.3399e-06	0.019236; ...
310	1.2726	472.98	0.0021142	295.73	298.43	1.3321	1.7212	2.8739	640.30	0.15178	8.5735e-05	0.088324	0.0055707	27.882	0.035865	566.42	612.07	2.3438	1.6531	2.1942	210.16	20.354	8.7494e-06	0.020863; ...
320	1.5992	454.98	0.0021979	324.18	327.70	1.4228	1.7654	3.0338	571.35	0.33451	7.6800e-05	0.083951	0.0044313	35.748	0.027973	575.25	619.98	2.3362	1.7180	2.3943	204.26	19.729	9.2300e-06	0.022776; ...
330	1.9835	434.92	0.0022993	353.94	358.50	1.5149	1.8134	3.2526	500.17	0.59757	6.8303e-05	0.079718	0.0033463	45.780	0.021844	582.96	626.29	2.3264	1.7922	2.6860	196.80	19.344	9.8173e-06	0.025121; ...
340	2.4320	411.87	0.0024280	385.45	391.35	1.6098	1.8670	3.5872	425.33	1.0100	6.0014e-05	0.075581	0.0023271	58.924	0.016971	588.97	630.25	2.3124	1.8839	3.1665	187.44	19.120	1.0574e-05	0.028187; ...
350	2.9527	383.88	0.0026050	419.52	427.21	1.7099	1.9335	4.2077	344.23	1.7479	5.1592e-05	0.071489	0.0013917	77.072	0.012975	592.13	630.44	2.2906	1.9998	4.1109	175.79	18.970	1.1635e-05	0.032683];

% isobaric properties - p = 1 Pa (for ideal gas cp)
% Data retrieved 2014-09-30, http://webbook.nist.gov/chemistry/fluid
%
%Temperature (K)	Pressure (MPa)	Density (kg/m3)	Volume (m3/kg)	Internal Energy (kJ/kg)	Enthalpy (kJ/kg)	Entropy (J/g*K)	Cv (J/g*K)	Cp (J/g*K)	Sound Spd. (m/s)	Joule-Thomson (K/MPa)	Viscosity (Pa*s)	Therm. Cond. (W/m*K)	Phase
p1Pa = ...
[120	1e-6	4.4196e-05	22627.	375.96	398.59	3.8699	0.81971	1.0083	166.83	233.45	3.4667e-06	0.0034998; ... vapor
130	1e-6	4.0796e-05	24512.	384.33	408.84	3.9519	0.85334	1.0419	173.00	181.75	3.7168e-06	0.0040849; ... vapor
140	1e-6	3.7882e-05	26398.	393.02	419.42	4.0303	0.88592	1.0745	178.93	144.77	3.9704e-06	0.0046992; ... vapor
150	1e-6	3.5357e-05	28283.	402.04	430.33	4.1055	0.91792	1.1065	184.64	117.55	4.2268e-06	0.0053426; ... vapor
160	1e-6	3.3147e-05	30169.	411.38	441.55	4.1780	0.94980	1.1384	190.15	97.010	4.4856e-06	0.0060152; ... vapor
170	1e-6	3.1197e-05	32054.	421.04	453.10	4.2479	0.98194	1.1705	195.47	81.164	4.7465e-06	0.0067169; ... vapor
180	1e-6	2.9464e-05	33940.	431.02	464.96	4.3158	1.0147	1.2032	200.62	68.710	5.0090e-06	0.0074478; ... vapor
190	1e-6	2.7913e-05	35826.	441.34	477.16	4.3817	1.0482	1.2368	205.60	58.761	5.2729e-06	0.0082078; ... vapor
200	1e-6	2.6517e-05	37711.	451.99	489.70	4.4460	1.0828	1.2714	210.42	50.698	5.5379e-06	0.0089969; ... vapor
210	1e-6	2.5255e-05	39597.	463.00	502.59	4.5089	1.1185	1.3070	215.11	44.082	5.8037e-06	0.0098153; ... vapor
220	1e-6	2.4107e-05	41482.	474.37	515.85	4.5706	1.1554	1.3439	219.66	38.595	6.0702e-06	0.010663; ... vapor
230	1e-6	2.3059e-05	43368.	486.11	529.48	4.6311	1.1934	1.3820	224.10	34.000	6.3371e-06	0.011539; ... vapor
240	1e-6	2.2098e-05	45253.	498.24	543.49	4.6908	1.2327	1.4212	228.42	30.119	6.6043e-06	0.012445; ... vapor
250	1e-6	2.1214e-05	47139.	510.77	557.90	4.7496	1.2731	1.4616	232.64	26.816	6.8715e-06	0.013380; ... vapor
260	1e-6	2.0398e-05	49024.	523.70	572.73	4.8077	1.3145	1.5030	236.76	23.988	7.1387e-06	0.014344; ... vapor
270	1e-6	1.9643e-05	50910.	537.06	587.97	4.8653	1.3568	1.5454	240.80	21.549	7.4056e-06	0.015337; ... vapor
280	1e-6	1.8941e-05	52796.	550.84	603.64	4.9222	1.4001	1.5886	244.76	19.436	7.6722e-06	0.016360; ... vapor
290	1e-6	1.8288e-05	54681.	565.06	619.74	4.9787	1.4440	1.6325	248.64	17.596	7.9382e-06	0.017411; ... vapor
300	1e-6	1.7678e-05	56567.	579.72	636.29	5.0348	1.4885	1.6771	252.45	15.985	8.2038e-06	0.018492; ... vapor
310	1e-6	1.7108e-05	58452.	594.83	653.29	5.0906	1.5336	1.7221	256.20	14.571	8.4686e-06	0.019602; ... vapor
320	1e-6	1.6573e-05	60338.	610.40	670.73	5.1460	1.5790	1.7675	259.89	13.325	8.7326e-06	0.020741; ... vapor
330	1e-6	1.6071e-05	62223.	626.41	688.64	5.2010	1.6246	1.8132	263.52	12.220	8.9957e-06	0.021909; ... vapor
340	1e-6	1.5598e-05	64109.	642.89	707.00	5.2559	1.6705	1.8590	267.11	11.238	9.2579e-06	0.023106; ... vapor
350	1e-6	1.5153e-05	65994.	659.82	725.82	5.3104	1.7163	1.9049	270.64	10.362	9.5190e-06	0.024333]; %	vapor

% Saturation pressure - psat < 80 Pa for T < 150
N.ps.Tfirst = 140;
N.v.Tlast = 330;
N.hvap.Tlast = N.v.Tlast;
N.jt.Tlast = N.v.Tlast;
N.cpg.Tlast = 310;
N.nug.Tlast = N.v.Tlast;

otherwise
error('Substance %s not known', name1);
end

s = substance(name1);
p = substance(name2);
allT = satdata(:,1);
len = size(allT,1);
first = 1;
last = len;

setfirstlast('ps');
nonarrayplot(s.ps,p.ps,'Saturation pressure',first,last,...
  allT,'T [K]',satdata(:,2)*1e6);

% Liquid density
setfirstlast('rho');
nonarrayplot(s.rho,p.rho,'Density of the liquid',first,last,...
  allT,'T [K]',satdata(:,3));

% enthalpy of vaporization
setfirstlast('hvap');
nonarrayplot(s.hvap,p.hvap,'Enthalpy of vaporization',first,last,...
  allT,'T [K]',(satdata(:,18)-satdata(:,6))*1e3);

% Specific volume of the vapor
setfirstlast('v');
pplot(s.v,p.v,'Specific volume of the vapor',first,last,...
  allT,satdata(:,2)*1e6,'T [K]',satdata(:,16));

% Joule-Thomson coefficient
setfirstlast('jt');
pplot(s.jt,p.jt,'Joule-Thomson coefficient',first,last,...
  allT,satdata(:,2)*1e6,'T [K]',satdata(:,23)*1e-6);

% Specific isobaric heat capacity of the vapor
setfirstlast('cpg');
pplot(s.cpg,p.cpg,'Specific isobaric heat capacity of the vapor',first,last,...
  allT,satdata(:,2)*1e6,'T [K]',satdata(:,21)*1e3);

% Specific isobaric heat capacity in the ideal gas state
setfirstlast('reset');
pplot(s.cpg,p.cpg,'Isobaric heat capacity, ideal gas state',1,size(p1Pa,1),...
  p1Pa(:,1),p1Pa(:,2)*1e6,'T [K]',p1Pa(:,9)*1e3);

% Specific isobaric heat capacity of the liquid
arrayplot(s.cpl,p.cpl,'Specific isobaric heat capacity of the liquid',first,last,...
  allT,'T [K]',satdata(:,9)*1e3);

% Surface tension
setfirstlast('sigma');
arrayplot(s.sigma,p.sigma,'Surface tension',first,last,allT,'T [K]',satdata(:,14));

% Dynamic viscosity of the liquid
setfirstlast('reset');
arrayplot(s.mul,p.mul,'Dynamic viscosity of the liquid',first,last,...
  allT,'T [K]',satdata(:,12));

% Dynamic viscosity of the vapor
arrayplot(s.mug,p.mug,'Dynamic viscosity of the vapor',first,last,...
  allT,'T [K]',satdata(:,24));

% Thermal conductivity of the vapor
setfirstlast('kg');
arrayplot(s.kg,p.kg,'Thermal conductivity of the vapor',first,last,...
  allT,'T [K]',satdata(:,25));

% Thermal conductivity of the liquid
setfirstlast('reset');
arrayplot(s.kl,p.kl,'Thermal conductivity of the liquid',first,last,...
  allT,'T [K]',satdata(:,13));

% Kinematic viscosity of the vapor
setfirstlast('nug');
pplot(s.nug,p.nug,'Kinematic viscosity of the vapor',first,last,...
  allT,satdata(:,2)*1e6,'T [K]',satdata(:,24).*satdata(:,16));

% Kinematic viscosity of the liquid
setfirstlast('reset');
nonarrayplot(s.nul,p.nul,'Kinematic viscosity of the liquid',first,last,...
  allT,'T [K]',satdata(:,12).*satdata(:,4));

function setfirstlast(var) %---- nested function --------
if isfield(N,var)
  if isfield(N.(var),'Tfirst')
    first = find(allT>N.(var).Tfirst,1);
  else
    first = 1;
  end
  if isfield(N.(var),'Tlast')
    last = find(allT>N.(var).Tlast,1) - 1 ;
  else
    last = len;
  end
else
  first = 1;
  last = len;
end
end %---------------------- end nested function ---------


function nonarrayplot(prop,qroq,description,first,last,X,xname,datanist) %------

global VERBOSE;
len = last - first + 1;
datam = zeros(len,1);
datan = datam;
for i = 1:len
  datam(i,1) = prop(X(first-1+i,1));
  datan(i,1) = qroq(X(first-1+i,1));
end

figure('Name',description,'NumberTitle','Off');
plot(X(first:last),100*(datam./datanist(first:last)-1),'ko',...
     X(first:last),100*(datan./datanist(first:last)-1),'kx');
legend('propane','propaneperry');
title(description);
ylabel('100*(matlab - nist)/nist');
xlabel(xname);

if VERBOSE
  printtable(description,X(first:last),datanist(first:last),datam);
end

end %---------------------------------------------------------------------------


function arrayplot(prop,qroq,description,first,last,X,xname,datanist) %---------

global VERBOSE;
datam = prop(X(first:last));
datan = qroq(X(first:last));
figure('Name',description,'NumberTitle','Off');
plot(X(first:last),100*(datam./datanist(first:last)-1),'ko',...
     X(first:last),100*(datan./datanist(first:last)-1),'kx');
legend('propane','propaneperry');
title(description);
ylabel('100*(matlab - nist)/nist');
xlabel(xname);

if VERBOSE
  printtable(description,X(first:last),datanist(first:last),datam);
end

end %---------------------------------------------------------------------------


function pplot(prop,qroq,description,first,last,X,P,xname,datanist) %-----------

global VERBOSE;
datam = prop(X(first:last),P(first:last));
datan = qroq(X(first:last),P(first:last));
figure('Name',description,'NumberTitle','Off');
plot(X(first:last),100*(datam./datanist(first:last)-1),'ko',...
     X(first:last),100*(datan./datanist(first:last)-1),'kx');
legend(name1, name2);
title(description);
ylabel('100*(matlab - nist)/nist');
xlabel(xname);

if VERBOSE
 printtable(description,X(first:last),datanist(first:last),datam,P(first:last));
end

end %---------------------------------------------------------------------------


function pscalarplot(prop,~,description,first,last,X,P,xname,datanist) %-----

global VERBOSE;
len = last - first + 1;
datam = zeros(len,1);
datan = datam;
for i = 0:len-1
  datam(i+1,1) = prop(X(first+i,1),P(first+i,1));
end

figure('Name',description,'NumberTitle','Off');
plot(X(first:last),100*(datam./datanist(first:last)-1),'ko',...
     X(first:last),100*(datan./datanist(first:last)-1),'kx');
legend(name1, name2);
title(description);
ylabel('100*(matlab - nist)/nist');
xlabel(xname);

if VERBOSE
 printtable(description,X(first:last),datanist(first:last),datam,P(first:last));
end

end %---------------------------------------------------------------------------


end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end testsubstance %%%

function printtable(description,X,datanist,datam,P) %---------------------------

fprintf(['\n' description ' \n']);
if nargin == 4
  fprintf('T [K]\tnist\tmatlab\t100*(matlab-nist)/nist\n');
  fprintf('%3.0f\t%.4g\t%.4g\t%.3f\n', ...
    [X'; datanist'; datam'; 100*(datam./datanist-1)']);
elseif nargin == 5
  if P(1) < 1e3
    fprintf('T [K]\tp [Pa]\tnist\tmatlab\t100*(matlab-nist)/nist\n');
    fprintf('%3.0f\t%.0f\t%.4g\t%.4g\t%.3f\n', ...
      [X'; P'; datanist'; datam'; 100*(datam./datanist-1)']);
  else
    fprintf('T [K]\tp [bar]\tnist\tmatlab\t100*(matlab-nist)/nist\n');
    fprintf('%3.0f\t%.4g\t%.4g\t%.4g\t%.3f\n', ...
      [X'; P'/1e5; datanist'; datam'; 100*(datam./datanist-1)']);
  end
end

end %---------------------------------------------------------------------------
