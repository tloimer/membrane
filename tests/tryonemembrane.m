addpath('../program');

T1 = 273.15 + 23;
s = substance('r142b');
%s = substance('isobutane');
f = fmodel('parallel');

% membrane(pore diameter, epsilon, thermal cond., tname, tau, beta, thickness)
mem = membrane(25e-9, 0.15, 1.22, 'tube', 1, 9.0541, 1e-4);

% mstackstruct(theta, {{layer1 layer2} {second membrane}}, fmodel)
ms = mstackstruct(0, {{mem}}, f);

p1 = 2e5;
p2 = 1e5;

[~,madi] = mnumadiabat(T1, p1, p2, s, ms);
ms.printsolution(madi);
ms.plotsolution(madi);

fprintf("Mnumadiabat succesful");
[~,mh] = mnumheatflux(30, T1, p1, p2, s, ms);
ms.printsolution(mh);
ms.plotsolution(mh);
