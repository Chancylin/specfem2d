## Note:

The work on this branch is developed by Chuangxin Lin (University of Toronto) as a conceptual study for the hybrid method. 

Publication (preprint version) is avaliable on researchgate webpage:

https://www.researchgate.net/profile/Kai_Wang128/publication/336286893_High-frequency_seismic_wave_modelling_of_the_deep_Earth_based_on_hybrid_methods_and_spectral-element_simulations_a_conceptual_study/links/5db895d94585151435d1651c/High-frequency-seismic-wave-modelling-of-the-deep-Earth-based-on-hybrid-methods-and-spectral-element-simulations-a-conceptual-study.pdf

Citation:

Lin, C., Monteiller, V., Wang, K., Liu, T., Tong, P. and Liu, Q., 2019. High-frequency seismic wave modelling of the deep Earth based on hybrid methods and spectral-element simulations: a conceptual study. Geophysical Journal International, 219(3), pp.1948-1969.



Instructions on how to install and use SPECFEM2D are
available in the PDF manual located in directory doc/USER_MANUAL.

Main "historical" developers: Dimitri Komatitsch and Jeroen Tromp
  (there are currently many more!)

For a quick test, run the default example with these commands:

  ./configure FC=gfortran
  make all
  ./bin/xmeshfem2D
  ./bin/xspecfem2D

and check the output files in ./OUTPUT_FILES/

