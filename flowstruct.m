%FLOWSTRUCT Structure containing the solution.
%  FLOWSTRUCT includes three fields, FLOWSTRUCT.INFO, FLOWSTRUCT.SOL and
%  FLOWSTRUCT.FLOW which are also structures.
%  FLOWSTRUCT.INFO includes
%    .kap   permeability of the membrane [m2]
%    .m     mass flux [kg/m2s]
%    .L     thickness of the membrane [m]
%    .T0    given temperature in front of the membrane [K]
%    .p0    given pressure in front of the membrane [Pa]
%    .dp    given pressure difference [Pa]
%    .ph    string denoting the 2-phase flow model
%
%  FLOWSTRUCT.SOL includes
%    .len   length of the structure FLOW
%    .a1    vapor volume fraction behind the condensation front
%    .q1    heat flux behind the condensation front [J/m2s]
%    .T0    temperature in front of the membrane [K]
%    .Te    temperature downstream of the membrane [K]
%    .p0    pressure in front of the membrane [Pa]
%    .pe    pressure downstream of the membrane [Pa]
%    .de    position of the evaporation front within the membrane [m]
%    .df    position of the condensation front [m]
%             (negative, if liquid film is present)
%
%  FLOWSTRUCT.FLOW includes
%    .z     z-coordinate
%    .T     temperature
%    .p     pressure
%    .q     heat flux
%    .a     vapor volume fraction
%    .x     vapor mass fraction
%    .color plotting specification (e.g.: '+k:')
