# FeXIII-coronal-density 

#### **Main Aim:** 
Two small example notebooks that can be used to infer FE XIII coronal densities using a Chianti calculation look-up table. Target lines are the infrared Fe XIII 1074.6nm and 1079.8nm pair. 

#### Separate notebooks particularized for distinct data types/sources exist for:

1. NSO's DKIST Cryo-NIRSP or DL-NIRSP slit and imaging spectro-polarimeters
2. HAO's MLSO CoMP/uCoMP imaging full-disk polarimetric coronagraph

uCoMP/CoMP data can be downloaded from [MLSO](https://mlso.hao.ucar.edu/mlso_data_calendar.php?calinst=ucomp).

Cryo-NIRSP or DL-NIRSP data can be downloaded from the [DKIST Datacenter](https://dkist.data.nso.edu/).

Some sample data that can be downloaded and unarchived directly into your repository clone that enables running the example notebooks [can be found here](https://www.dropbox.com/scl/fi/b5flnx7jkyo2rb561cg8p/FeXIII-coronal-density_sampledata.tar.gz?rlkey=g7cawm4m1wzn3hkri04u1xkrw&dl=1). 

### **Notes:**
- This functionality is included in the CLEDB coronal inversion. This repository provides a standalone method for computing coronal densities only.
- Spatio-temporal matching observations in both Fe XIII infrared lines are required.  
- Examples offer a parallel cpu implementation making calculations rasonably fast.
- Options for processing either peak line emission or integrated line emission exist for spectroscopic instruments.
- Accurate header and pointing information from the data providers are required.
- Currently, a high-resolution look-up table is provided in this repository. The table is created using [CHIANTI V10.1](https://download.chiantidatabase.org/CHIANTI_10.1_database.tar.gz) through a [PyCELP](https://github.com/tschad/pycelp) implementation.
- A PyCELP look-up table generator notebook is also provided, but a user should not generate new look-up tables unless a specific need exists. The provided tables are flexible enough to suit most needs.
- A low-res older and deprecated look-up table created via a [SSWIDL](https://www.mssl.ucl.ac.uk/surf/sswdoc/solarsoft/ssw_setup.html) CHIANTI implementation exists, but this should only be used for cross-validating calculations, and not for production runs. 

#### **Contact:** Alin Paraschiv, NSO arparaschiv at nso edu

#### **Acknowledgements/Credits**

This implementation benefits and uses published results from:

1. [Paraschiv & Judge, SolPhys, 2022](https://ui.adsabs.harvard.edu/abs/2022SoPh..297...63P/abstract) 
2. [Schad & Dima, SolPhys, 2020](https://ui.adsabs.harvard.edu/abs/2020SoPh..295...98S/abstract)