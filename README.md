# 2D-PT: Two-Dimensional Parallel Tempering for Constrained Optimization

## Scope

This repository contains a MATLAB implementation of **2D-PT**, a novel extension of Parallel Tempering for constrained optimization. In addition to the temperature axis, 2D-PT introduces a new penalty dimension for constraints where replicas have increasing penalty strengths and exchange their configuration based on constraint satisfaction. For more details, please see our paper on arXiv: https://arxiv.org/pdf/2506.14781

The method is demonstrated on **sparsified Ising problems** where auxiliary copy spins are introduced to reduce the number of neighbors per spin, useful for practical hardware implementation with limited connectivity. To meet the constraints, copy spins should be positively coupled with a sufficient strength, at the risk of rigidifying the system and challenging the ground state search. Instead of fine-tuning the copy strength, 2D-PT explores multiple copy strengths (columns) and exploits feasible samples from the relaxed penalty strengths (lower columns) to speed up the constrained search.

## Content

- MATLAB implementation of 2D Parallel Tempering for sparsified graphs (`Source_code/dynamics_2D_PT_sparsified_graph.m`) with adaptive schedules (`Source_code/APT_2D_for_sparsified_graph.m`).
- 2D-PT example with sparsified Wishart instances (`Example_2DPT/main_2DPT_for_sparsified_Wishart.m`).
- Data and scripts to plot the manuscript figures (`Generate_figures/generate_figure.m`).

## Requirements

- Tested on MATLAB R2021 or later
- No external toolboxes required

## Contributing

Contributions to improve the code or extend its functionality to other constrained optimization problems are welcome. Please feel free to submit issues or pull requests.

## Acknowledgements

- Wishart instances were generated using Chook: https://github.com/dilinanp/chook
- Data is plotted in MATLAB using `shadedErrorBar`: https://github.com/raacampbell/shadedErrorBar

## Contact

Please feel free to contact Corentin Delacour (delacour@ucsb.edu) for any questions or suggestions.
