This is the repository containing all supplementary files for my [thesis](Thesis_ModelApproximation_REP.pdf) *Model Approximation in Transport Dynamics Motivated by Biological Phenomena*, 
written in fulfillment of the Bachelor of Science in Mathematics at Leiden University.

The folder [Figures](figures) contains a .drawio file that can be loaded into draw.io to edit the schematic figures presented in the thesis. Additionally it contains the folders [Snapshots_Nh_1000](Figures/Snapshots_Nh_1000), [Snapshots_Nh_2500](Figures/Snapshots_Nh_2500) and [Snapshots_Nh_5000](Figures/Snapshots_Nh_5000) that contain figures similar to those presented in Figure 7 in my thesis, with the number of hitchhikers as displayed in the directory name.

The .R-files contain code used in the simulations, with the parameter values defined in [definitions_VJ.R](definitions_VJ.R) and the numerical schemes in [definitions_numerics.R](definitions_numerics.R). The code in [calculate_D.R](calculate_D.R) is the implementation of Algorithm 1 presented in Section 7.4 of the thesis. The code in [solve_pde_u.R](solve_pde_u.R) and [comparison_unif2.R](comparison_unif2.R) contain the code used to compare the particle-based model with the deterministic PDE-model of the Hitchhiking phenomenon, as described by Algorithm 2 in the same section. 

For the simulation of the moments-based system described in Section 9, the file [moments_system.R](moments_system.R) was used.

The ensemble-MSD, was computed using [msd_experiment.R](msd_experiment.R), and the single particle MSD was simulated using [msd_single_particle.R](msd_single_particle.R). The two calculated MSDs were compared using [plot_msd.R](plot_msd.R), as shown in Figure 10 in Section 9.
