%FLOWSTRUCT Structure containing the solution.
%  FLOWSTRUCT includes four fields, FLOWSTRUCT.INFO, FLOWSTRUCT.SOL,
%  FLOWSTRUCT.FLOW and FLOWSTRUCT.LIN which are also structures.
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
%    .a3    vapor volume fraction behind the condensation front
%    .q3    heat flux behind the condensation front [W/m2]
%    .T0    temperature in front of the membrane [K]
%    .Te    temperature downstream of the membrane [K]
%    .p0    pressure in front of the membrane [Pa]
%    .pe    pressure downstream of the membrane [Pa]
%    .de    position of the evaporation front within the membrane [m]
%    .df    position of the condensation front [m]
%             (negative, if liquid film is present)
%
%  FLOWSTRUCT.FLOW includes
%    .z     z-coordinate [m]
%    .T     temperature [K]
%    .p     pressure [Pa]
%    .q     heat flux [W/m2]
%    .a     vapor volume fraction
%    .x     vapor mass fraction
%    .color plotting specification (e.g.: '+k:')
%
%  FLOWSTRUCT.LIN gives the solution according to the linear theory,
%    .m     mass flux [kg/m2s]
%    .Te    temperature at the evaporation front [K]
%    .T4    temperature at the front of the membrane [K]
%    .p6    pressure at the evaporation front [Pa]
%    .deL   relative position of the evaporation front [-]
%    .dfL   relative position of the film [-]
%    .a3    vapor volume fraction behind the condensation front
%    .x3    vapor mass fraction behind the condensation front
