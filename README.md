# FeXIII-coronal-density 

#### **Main Aim:** 

This repository provides a straightforward and standalone implementation for computing coronal densities.

Two small example notebooks that can be used to infer FE XIII coronal densities using a CHIANTI calculation look-up table. Target lines are the infrared Fe XIII 1074.7nm and 1079.8nm pair. 

#### **Two notebooks particularized for distinct data types/sources exist for:**

1. NSO's DKIST Cryo-NIRSP or DL-NIRSP slit and imaging spectro-polarimeters
2. HAO's MLSO CoMP/uCoMP imaging full-disk polarimetric coronagraph

#### **Sample data download:**

Some sample data that can be downloaded and unarchived directly into your repository clone will enable running the example notebooks [can be found here](https://drive.google.com/file/d/1I0vn7tJaMn4EU_RK3b1YxMYpPmiUZsue/view?usp=drive_link). 

uCoMP/CoMP data can be downloaded from [MLSO](https://mlso.hao.ucar.edu/mlso_data_calendar.php?calinst=ucomp).

Cryo-NIRSP or DL-NIRSP data can be downloaded from the [DKIST Datacenter](https://dkist.data.nso.edu/).

#### **Notes, Assumptions, and Caveats:**
- This implementation is a standalone coronal density calculations that are implemented and available through the [CLEDB coronal inversion](https://github.com/arparaschiv/solar-coronal-inversion).
- Spatio-temporal matching observations in both Fe XIII infrared lines are required.
- **The LOS integration assumption exists.** Users should be very careful when interpreting their coronal signals as they get integrated along the LOS. This approach is valid only in spaces where a dominant coherent structure exists. A superposition of structures along the LOS is not adequately represented in terms of density computed via this method.
- Examples offer a parallel/multithread cpu implementation, making calculations reasonably fast.
- Options for processing either peak line emission or integrated line emission exist for all instruments.
- Accurate header and pointing information from the data providers are required. **Currently, the DKIST Cryo-NIRSP .asfd metadata has known bugs. Please be mindful of pointing, and recreate a more feasible spatial coordinate matching.**
- Currently, a high-resolution look-up table is provided in this repository. The table is created using [CHIANTI V10.1](https://download.chiantidatabase.org/CHIANTI_10.1_database.tar.gz) through a [PyCELP](https://github.com/tschad/pycelp) implementation.
- A PyCELP generated look-up table generator notebook is also provided. Users should not generate new look-up tables unless a specific need exists. The provided tables are flexible enough to suit most needs.
- A lower parameter space resolution deprecated look-up table created via a [SSWIDL](https://www.mssl.ucl.ac.uk/surf/sswdoc/solarsoft/ssw_setup.html) CHIANTI implementation still exists, but this is deprecated and should only be used for cross-validating calculations, and not for production runs. 

#### **Contact:** Alin Paraschiv, NSO  -- arparaschiv at nso edu

#### **Acknowledgements/Credits**

This implementation benefits and uses published results from:

1. [Paraschiv & Judge, SolPhys, 2022](https://ui.adsabs.harvard.edu/abs/2022SoPh..297...63P/abstract) 
2. [Schad & Dima, SolPhys, 2020](https://ui.adsabs.harvard.edu/abs/2020SoPh..295...98S/abstract)
3. [Del Zanna, G., et. al, ApJ, 2021](https://ui.adsabs.harvard.edu/abs/2021ApJ...909...38D)
