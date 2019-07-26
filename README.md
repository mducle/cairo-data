# Inelastic neutron data on cairo-lattice materials

This repository contains the raw inelastic neutron scattering data files for experiments on the Cairo lattice materials Bi<sub>2</sub>Fe<sub>4</sub>O<sub>9</sub> and Bi<sub>4</sub>Fe<sub>5</sub>O<sub>13</sub>F.
It also has python files which will process the raw data files into a format readable by the standard data visualisation package, MSlice, used for these data. 
There are three incarnations of this program, [the original](http://mslice.isis.rl.ac.uk/Main_Page) which requires Matlab, an [IDL rewrite](https://www.ncnr.nist.gov/dave/download.html) within the DAVE program, and a [Python rewrite](https://mantidproject.github.io/mslice/) within the [Mantid](https://download.mantidproject.org/) program.
Mantid is required to run the reduction scripts, but pre-processed reduced data files are also given in this archive.
