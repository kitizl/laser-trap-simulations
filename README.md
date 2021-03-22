# Laser Trap Simulations

A repository containing simulations and statistics to study nanosphere dynamics in laser traps

---


# Download

Download the repository using `git clone https://github.com/kitizl/laser-trap-simulations`. You will need Julia to run the .jl files and Python to run the .py files. This is meant for scientific computation, so it is assumed that you have a sufficiently up to date scientific package distribution, such as Anaconda.

# Usage
The structure of this repository will soon change and therefore, the usage documentation will be updated then.

## Single-well

The key interaction occurs through `driver.py`, dedicated to running the simulations and statistics based on a specific setup. Through command line arguments, you are able to set the temperature, pressure, and other parameters for the simulations. For more in depth fine tuning of the parameters, such as the mass and volume of the nanosphere, and the specifics of the stiffness of the trap, they still need to be changed, either in `simulation.py`, where the key physics lies or in `statistics.py` where the analysis of the data is stored.

```bash
Usage: python driver.py [-h] -t MAX_TIME -p PRESSURE -T TEMPERATURE -s SAVING_FREQ -N NUMTRIALS -r RESOLUTION [-l LABEL]

optional arguments:
  -h, --help            show this help message and exit
  -l LABEL, --label LABEL
                        label for the experiment
mandatory arguments:
  -t MAX_TIME, --max_time MAX_TIME
                        maximum time for simulation (in s)
  -p PRESSURE, --pressure PRESSURE
                        pressure (in mbar)
  -T TEMPERATURE, --temperature TEMPERATURE
                        temperature (in K)
  -s SAVING_FREQ, --saving_freq SAVING_FREQ
                        saving frequency (number of steps per save)
  -N NUMTRIALS, --numTrials NUMTRIALS
                        number of trials that need to be run
  -r RESOLUTION, --resolution RESOLUTION
                        resolution for the histogram
```

The data will be stored in a folder prefixed `data`, followed by either the label you have provided, or a timestamp corresponding to when the run began. The plots produced using `statistics.py` as called from `driver.py` will be stored in a folder prefixed `plot`, followed by either the label you have provided or a timestamp correspodning to when the run began.

## Double-well

This folder is not under active development, and might be removed in a future update when the single-well code can be (and as it stands, it in principle can) utilized to model a double well potential.