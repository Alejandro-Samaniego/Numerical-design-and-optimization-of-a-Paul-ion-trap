# Simulating a Paul Trap
## Project

This project was carried out as the final project on the course PEF2 and it was developed by Janot Vilaró(https://github.com/janotvilaro) and me.

## Overview
We expose the theoretical concepts and numerical tools required to simulate a Paul trap. We demonstrate that it is possible to confine ions by designing our own RF trap and we focus our study on $^9\text{Be}^+$ ions which are of great interest for their applications in quantum information processing and quantum computing [Kielpinski et al., 2002; Monroe et al., 1995].

Additionally, we perform a deeper study for the confinement of two ions: we highlight the relevance of the relative masses ($m_{\text{rel}} = \frac{m_2}{m_1}$) and charges ($q_{\text{rel}} = \frac{q_2}{q_1}$), together with a rough stability test of the frequency-voltage relation that assures optimal confinement.

## Ion Confinement
<p align="center">
  <img src="Figures/2 ions.jpg" alt="Figure 1" width="400"/>
  <img src="Figures/7 ions.jpg" alt="Figure 2" width="400"/>
</p>

<p align="center">
  <gif src="Figures/Paul Trap gif.gif alt="Figure 3" width="400"/>
</p>

## Comments on use

Make sure to download all of the files of the "Code" section and save them on your computer in a single folder. 

Open the file named PaulTrap. The work is divided in sections, which target increasingly harder problemes. Consider running them ordely as the value of some variables may be needed for later sections. The end of the document performs a more detailed analysis of convergence for the 2-ion scenario.

Note that most of the .m files are functions used to generate the mesh on the hyperbolic surface, to integrate required equations or to integrate the ions trajectories, consequnetly you won't need to open nor execute such files.

## Convergence analysis

![Frequency vs Potential](Figures/stabilityfinale.jpg)

## References

- Kielpinski, D., Monroe, C., & Wineland, D. J. (2002). Architecture for a large-scale ion-trap quantum computer. *Nature*, 417(6890), 709-711.
- Monroe, C., Meekhof, D. M., King, B. E., Jefferts, S. R., Itano, W. M., Wineland, D. J., & Gould, P. (1995). Demonstration of a Fundamental Quantum Logic Gate. *Physical Review Letters*, 75(25), 4714-4717.
