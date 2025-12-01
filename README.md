# SPPARKS_TMA-H2O
SPPARKS for TMA&amp;H2O ALD
KMC Model for Al₂O₃ ALD using Customized SPPARKS
This repository contains a customized SPPARKS implementation developed for simulating trimethylaluminum (TMA)–based atomic layer deposition of Al₂O₃ using a lattice kinetic Monte Carlo model. The code modifications and example files allow readers to reproduce the simulations presented in our manuscript.

Quick Start Guide  
----------------
Copy the customized source files  
Copy all files in the src/ directory of this repository into the src/ directory of your local SPPARKS installation.  

Compile SPPARKS  
Build SPPARKS as usual to generate the spk executable.  
(If you are unfamiliar with compilation, please refer to the official SPPARKS documentation：https://spparks.github.io)  

Run the provided example  
Navigate to the example/ folder and run the simulation using your compiled spk executable.


Notes  
----------------
The input script in.ald contains reaction energies corresponding to 393 K.

The provided structure file data.ald represents a 10×10 surface rather than the 20×20 system used in the manuscript, in order to reduce computational cost and file size.

An additional instruction.pdf is included to guide users on how to manually modify the in.ald file for different temperatures or surface configurations.

For more detailed explanations of syntax, commands, and KMC algorithms, please consult the official SPPARKS manual.
