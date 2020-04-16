# Some fancy LAr simulation for checking efficiency of LAr veto in L200, scanning attenuation lengths while still looking a bit into WLSR dimensions but never touching LAr veto dimensions. That we will do for L1000 ;)

TODO: 1. Change title.

My Idea for naming/simulating convention:

0. Simulation Campaign: set of several simulation programs running in parallel (or serial) and scanning individual voxels of the geometry for one fixed attenuation/WLSR size/LAr veto (the latter we swore to never touch)
1. RunSeries: All which is covered in one single program execution. I suspect to pack several voxels into one RunSeries to minimize overhead to start a large number of single jobs.
2. Run: =G4Run, i.e. one execution of /run/beamOn with a fixed nr of events to be simulated. ParticleGun setting (--> position, i.e. voxel) is fixed (save maybe for direction change or potential distribution within one voxel) as well as all the rest. Sure we don't touch LAr veto here either.
3. Event: =G4Event, i.e. all stuff simulated after one single particle gun shot. This is one single 128 nm photon in our case.
4. Hit: A single hit of the LAr veto instrumenation. By construction: an Event leads to 0 or 1 hit. Want to save position by giving a volumeID. cumulative hits in a single run are added up and stored at a voxel position.



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
