Delta 3D printer calibration utility in Python.
It is based on David Crocker's (DC42) javascript code which, in turn, implements the same algorithm found in his RepRapFirmware fork for deltas.

This is a command line utility aimed to:
1- calibrate a delta which have no z-probe
2- be scripted (e.g.: integrated in OctoPrint) to automate the delta calibration by offloading the process from the controller to a PC or Raspberry.