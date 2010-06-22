% Maximaler Fehler mit fehlerhafter Verwendung von drho in mnum, flow12:
% max(abs(mneu/mcalc-1)) = 0.0023

load ../matlab/results1006;
mneu=mcalc;
load ../matlab/results0906;
isvapor = ~strcmp(substancename,'nitrogen');

max(abs(mneu(~isvapor)./mcalc(~isvapor)-1))

max(abs(mneu(isvapor)./mcalc(isvapor)-1))
