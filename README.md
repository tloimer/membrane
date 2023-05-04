Matlab code to compute the flow of several substances, e.g., nitrogen, ethanol,
butane, isobutane, through homogeneous membranes, membranes with several,
different layers, and through stacks of several, individual membranes. The
membranes in a membrane stack may each consist of different layers.

The state of the substances upstream of the membrane, or membrane stack, must be
gaseous, and can be close to or at saturation. For the flow through the
membrane, condensation and evaporation due to capillary condensation or due to
the Joule-Thomson effect is taken into account. The energy balance and heat
transfer from or to the downstream side of a membrane is taken into account.

Make the program subfolder be in the path of matlab, e.g., by starting matlab
there. Type 'help program' to have an overview of the functions provided in
there.

The code is licensed under the
Creative Commons Attribution 4.0 International License.

For condensation due to the Joule-Thomson effect, see  
  W. Schneider. Vapor flow through a porous membrane — a throttling process
  with condensation and evaporation, _Acta Mechanica 47_, 15–25 (1983).
  [doi:10.1007/BF01176497](https://doi.org/10.1007/BF01176497).  
and  
  T. Loimer. Linearized description of the non-isothermal flow of a saturated
  vapor through a micro-porous membrane, _J. Membr. Sci. 301_, 107–117 (2007).
  [doi:10.1016/j.memsci.2007.06.005](https://doi.org/10.1016/j.memsci.2007.06.005).
  

Publications, for which this code was used, include

Thomas Loimer, Stepan K. Podgolin, Javad Sodagar-Abardeh,
Dmitrii I. Petukhov, Andrei A. Eliseev.
Influence of heat transfer and wetting angle on condensable fluid flow through
nanoporous anodic alumina membranes, _Phys. Chem. Chem. Phys. 25_, 3240–3250 (2023).
[doi:10.1039/d2cp04577j](https://doi.org/10.1039/d2cp04577j).

Katerina Setnickova, Roman Petrickovic, Petr Uchytil, Thomas Loimer.
Experimental and numerical study of the flux of isobutane vapors near
saturation through multi-layered ceramic membranes, _Sep. Purif.
Techn. 306_, 122604 (2023).
[doi:10.1016/j.seppur.2022.122604](https://doi.org/10.1016/j.seppur.2022.122604).

T. Loimer. The curvature of an evaporating meniscus in a pressure
driven flow through cylindrical pores, _Proc. Appl. Math. Mech. 19_,
e201900114 (2019).
[doi:10.1002/pamm.201900114](https://doi.org/10.1002/pamm.201900114).

T. Loimer, K. Setnickova, P. Uchytil. Consideration of the
Joule-Thomson effect for the transport of vapor through anodic alumina membranes
under conditions of capillary condensation, _Sep. Purif. Techn. 215_, 548–556 (2019).
[doi:10.1016/j.seppur.2019.01.051](https://doi.org/10.1016/j.seppur.2019.01.051).

P. Uchytil, J. Reznickova, K. Setnickova, T. Loimer. Comparison
of the flow of permanent and condensable gases through an asymmetric porous
membrane, _Chemie Ingenieur Technik 88_, 1779–1787 (2016).
[doi:10.1002/cite.201600047](https://doi.org/10.1002/cite.201600047).

T. Loimer, P. Uchytil. Influence of the flow direction on the mass
transport of vapors through membranes consisting of several layers. _Exp.
Thermal Fluid Sci. 67_, 2–5 (2015).
[doi:10.1016/j.expthermflusci.2014.12.012](https://doi.org/10.1016/j.expthermflusci.2014.12.012).

P. Uchytil, T. Loimer. Large mass flux differences for opposite flow
directions of a condensable gas through an asymmetric porous membrane. _J.
Membr. Sci. 470_, 451–457 (2014).
[doi:10.1016/j.memsci.2014.07.055](https://doi.org/10.1016/j.memsci.2014.07.055).

T. Loimer. The thermodynamic states of a fluid in a Joule-Thomson
process involving phase changes at interfaces with large curvature. In M.
Pilotelli, G. P. Beretta (Eds.). Proceedings of the 12th Joint European
Thermodynamics Conference, Brescia, Italy, July 1–5, 2013, 537–541 (Snoopy,
2013).
