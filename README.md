# Counterion-controlled phase equilibria in a charge-regulated polymer solution

This repository contains the JULIA codes used to study the phase equilibria of charge-regulated polymers in an acidic solution -- see discussion in the pre-print *Counterion-controlled phase equilibria in a charge-regulated polymer solution* available on the arXiv at the link https://arxiv.org/abs/2307.03706. 

* `Homogeneous_equilibria.jl` Jupyter notebook that computes the homogeneous equilibria of the mixture for the charge-regulated (CR) polymer model. The code generates the data used to generate *Figure 4* in the <a target="_blank" href="https://arxiv.org/abs/2307.03706">paper</a>.
* `Compute_phase_diagram_FC.jl` Jupyter notebook that computes the phase diagrams of the mixture for the fixed charge (FC) polymer model. The code generates the data used to generate *Figure 5* in the  <a target="_blank" href="https://arxiv.org/abs/2307.03706">paper</a>.
* `Compute_phase_diagram_CR.jl` Jupyter notebook that computes the phase diagrams of the mixture for the charge-regulated (CR) polymer model. The code generates the data used to generate *Figure 6* in the  <a target="_blank" href="https://arxiv.org/abs/2307.03706">paper</a>.

**NOTE:** The code requires downloading the latest version of the Julia package  <a target="_blank" href="https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/"> BifurcationKit.jl</a>
