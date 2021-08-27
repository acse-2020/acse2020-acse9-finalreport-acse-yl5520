# acse2020-acse9-finalreport-acse-yl5520

This repository contains the codes developed by Yuanshen Liao in his MSc. Individual Research Project.

## Usage

All required third-party modules to use the developed code are listed in `requirements.txt`

The `numpy.f2py`  is used to compile the FORTRAN library used in the developed code through issuing command: `f2py -c src/diffusion3d_virus_town_large.f90 -m diffusion3d_virus_town_large` under the root directory of this repository.

Two folders `data` and `images` are needed to create in advanced.

Command `python main.py --output save0` starts the simulation implemented by the developed code with the configuration files under `conf` directory, and output a file `save0.pkl` in `data/`

Command `python main.py --plot data/save0.pkl` plots all simulation results.

`python main.py -h` shows more detail usage of the developed code



## Configuration files

All six `JSON` files under `conf/` are loaded in `main.py` by default. Users can specify their own configuration file through the corresponding flag, see `python main.py -h` for more detail.

All keys in these six `JSON` files should remain the same in any alternative configuration files. Most of these keys are self explanatory, see below for explaination of some keys.

`distribution.json` defines mean value and standard derivation of normal distribution as a list of length 2. The first value is interpreted as mean value, the second value is interpreted as standard derivation.

`hosts-conf-default.json` defines the number of people of each attribute. See `ACSE_FinalReport_Virus_Yuanshen_Liao.pdf` for explanation of attributes.

`schedule-default.json` defines the time for people being inside or outside of the classroom. `"type": 1` indicates people are inside the classroom during the given period of time and are outside the classroom during the non-given period of time, `"type": 0` is vice versa. `"start"` indicates the starting time of the simulation, `"details"` indicates the periods of time. `"start"` and `"details"` must be in the same format. Keys in `"details"` must be the first three characters of a day in a week.



## Others

`scripts/` contains some python scritps used in development for analysing or plotting simulation output, it is not relevant to the simulation process.
