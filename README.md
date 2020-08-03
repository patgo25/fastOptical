# A fast optical simulation to track photons in the L200 geometry

## Goals

- Sweeping optical parameters/geometries in minutes
- Utilizing symmetry assumptions
- Implementing fancy tricks to speed up simulation

-> A grey area between the black and white of super precise Monte Carlo MaGe (slow) and back of the envelope calculations. Hence not a replacement for MaGe. First find best fitting parameters with this simulation and then run MaGe once

## Symmetry assumptions

A rotational symmetry of 360/28 deg is assumed (half of a detector string). 4 scan modes (from fastest to slowest runtime):
- 1D scan: along the x-axis at a defined z-position
- 2D scan: Scans area to obtain a full map (360/28 deg of the cake) at a defined z-position
- XZ scan: Scans the XZ plane at a defined y-position
- 3D scan: Scans area to obtain a full map (360/28 deg of the cake) and a defined z-range

## Fancy tricks

- One of the main reasons why optical simulations are so slow in GEANT4 are the optical border processes. Hence here the amount of opical borders is minimized and at the ctrical places (e.g. LAr volume <-> fiber volume) a costum physics process is implemented.

- Another CPUh black hole is the propagation of the photons in the fiber. This is in this simulation completely avoided by doing an analytical propagation with parameters based on measurements

Again one has to emphasize that this simulation approach can not replace a full Monte Carlo, since all this tricks introduce slight errors from second order processes (e.g. photon leaves fiber during the propagation and couples into another fiber and gets detected. While the analytical model accounts for photons leaving a fiber it does not account for these photons beeing able to couple back into another fiber)

## Results
- Speed: 20000 - 60000 events/s (Intel(R) Xeon(R) CPU E5-2640 v3 @ 2.60GHz) depending on the choosen optical parameters (mainly attenuation length of the LAr). For compharsion MaGe: 500 - 1500 events/s (Intel(R) Xeon(R) CPU E5-2643 v3 @ 3.40GHz)
- Accuracy: Simulation results agree with the MaGe results

## Future

Harness the power of GPUs which are thanks to the video game industry optimiced for ray tracing?

----

And please use cmake. the makefile is unmaintained as of now

 based on g4simple, see here:

# g4simple
Perhaps the simplest fully-featured G4 application.

Based on one-file simulation by Jason Detwiler.

Installation:
Compile geant4 with GDML support (and optionally HDF5 support), then do:

```source (g4install_path)/share/(g4version)/geant4make/geant4make.sh```

(or `source (g4install_path)/share/(g4version)/geant4make/geant4make.csh`)

```make```

Physics List: uses Geant4's 
[named physics lists](https://geant4.web.cern.ch/node/155), 
set them using macro commands (see example run.mac)

Generator: uses Geant4's 
[GPS](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/GettingStarted/generalParticleSource.html). 
Set it up using macro commands (see example run.mac).

Geometry: uses 
[GDML](http://lcgapp.cern.ch/project/simu/framework/GDML/doc/GDMLmanual.pdf) 
(see example run.mac). Can use materials from Geant4's
[NIST Material Database](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/Appendix/materialNames.html) (note: the GDML parser will complain that the materials have not been defined, but Geant4 will still run without error).
Also supports Geant4's [text file geometry scheme](http://geant4.cern.ch/files/geant4/collaboration/working_groups/geometry/docs/textgeom/textgeom.pdf).

Output: uses Geant4's analysis manager (root, hdf5, xml, csv), with several
configurable options for output format, sensitive volumes (including [regex-based pattern matching / replacement](http://www.cplusplus.com/reference/regex/ECMAScript)), etc. (see example
run.mac). Records event/track/step numbers, 
[PIDs](http://pdg.lbl.gov/2018/reviews/rpp2018-rev-monte-carlo-numbering.pdf) 
(see also the python package [particle](https://pypi.org/project/Particle/)),
positions, energies, etc.

Visualization: uses avaialable options in your G4 build (see example vis.mac).

Postprocessing: you will want to postprocess the output to apply the detector
response. See example code that runs on the output of run.mac.


See similar project by Jing Liu at https://github.com/jintonic/gears .
