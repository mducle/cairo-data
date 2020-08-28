# Inelastic neutron data on cairo-lattice materials

This repository contains the raw inelastic neutron scattering data files for experiments on the Cairo lattice materials Bi<sub>2</sub>Fe<sub>4</sub>O<sub>9</sub> and Bi<sub>4</sub>Fe<sub>5</sub>O<sub>13</sub>F.
It also has python files which will process the raw data files into a format readable by the standard data visualisation package, MSlice, used for these data. 
There are three incarnations of this program, [the original](http://mslice.isis.rl.ac.uk/Main_Page) which requires Matlab, an [IDL rewrite](https://www.ncnr.nist.gov/dave/download.html) within the DAVE program, and a [Python rewrite](https://mantidproject.github.io/mslice/) within the [Mantid](https://download.mantidproject.org/) program.
Mantid is required to run the reduction scripts, but pre-processed reduced data files are also given in this archive.

In addition, the `calculations` folder contains Matlab data analysis scripts for fitting linear spin wave models using the [SpinW](https://spinw.org/) program to the data.
Also included are input files for McPhase calculation of the single-ion anisotropy using a point charge model and additional mean-field modelling of the high field magnetisation.

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work by Manh Duc Le is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
