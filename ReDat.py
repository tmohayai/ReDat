#!/usr/bin/env python

"""
Author: Tanaz A. Mohayai
Copyright (C) 2016-present by Tanaz A. Mohayai. All rights reserved.

This module is tailored for use with a Muon Ionization Cooling Experiment (MICE) data file in the ROOT format. It reads this
ROOT file which has a MICE-specific data structure and converts it into a for009 formatted file for input into ecalc9f for 
further analysis. 

This script, along with others serve the Muon Ionization Cooling Experiment, MICE and are part of the author's PhD thesis 
work on MICE [1, 2, 3, 4, 5]. The function parameters defined in this module along with their descriptions are as following: 

root_file_name: name of the input ROOT file.

tof_min, tof_max: user-defined range for the time-of-flight of muons. 

plane_number: MICE trackers consist of 5 scintillating fiber stations, each with 3 planes. User can specify a particular 
tracker plane for emittance or density calculation. 

usp_min, usp_max: user-defined range for the tracker reconstructed muon momenta. Together with TOF time-of-flight, they act 
as PID cuts.  

Note: only relevant columns of the for009 file are populated. MICE data does not have information regarding some of the columns 
(polarization, weight, etc) in the for009 file. 

Please note that this module, along with others which are currently in author's posession are planned to be merged into 
the official MICE Analysis User Software, MAUS package [8].

[1] Scipy's "gaussian_kde()" module by R. Kern, http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html.
[1] T. A. Mohayai et al., "Novel Application of Non-parametric Density Estimation Technique in Muon Ionization Cooling Experiment," APS-DPF'17 Proceedings.
[2] T. A. Mohayai, "Novel Application of Kernel Density Estimation in MICE," MICE-Note-506.
[3] T. A. Mohayai et al., "Novel Implementation of Non-parametric Density Estimation in MICE," IPAC'17 Proceedings.
[4] T. A. Mohayai et al., "Simulated Measurements of Beam Cooling in Muon Ionization Cooling Experiment," NA-PAC'16 Proceedings.
[5] T. A. Mohayai et al., "Simulated Measurements of Beam Cooling in Muon Ionization Cooling Experiment," NA-PAC'16 Proceedings.
[6] T. Roberts, "G4beamline User's Guide", Muons, Inc (2013), http://www.muonsinc.com. 
[7] http://www.cap.bnl.gov/ICOOL/fernow/readme.html
[8] C. D. Tunnell, C. T. Rogers, "MAUS: MICE Analysis User Software", IPAC (2011).
"""

import os
import subprocess
# basic PyROOT definitions
import ROOT 
# definitions of MAUS data structure for PyROOT
import libMausCpp 
import numpy as np
from math import *
import itertools

def redat(root_file_name, tof_min, tof_max, plane_number, usp_min, usp_max):
	
	# loop variable used for proper event number assignments in upstream and downstream trackers.
	p_0=0
	
	usfile=open('for009_US.dat','w')
	# the column titles are selected for illustration-purposes only; only the MICE data relevant for transverse phase-space calculations have titles.
	usfile.write('#  evt region x y z px py pz station_number\n')

	dsfile=open('for009_DS.dat','w')
	dsfile.write('#  evt region x y z px py pz station_number\n')

	# read the root file.
	root_file = ROOT.TFile(root_file_name,"R")
	# MICE-specific PyRoot commans; for further information, please refer to MAUS user guide.
	data = ROOT.MAUS.Data() # pylint: disable = E1101
	tree = root_file.Get("Spill")
	tree.SetBranchAddress("data", data)
	# iterate over the number of spills
	for i in range(tree.GetEntries()):
		tree.GetEntry(i)
		spill = data.GetSpill()
		# data quality cut: only physical events are allowed
		if spill.GetDaqEventType() == "physics_event":
			
			for j in range(spill.GetReconEvents().size()):
				p_0 = p_0+1
				recon_event=spill.GetReconEvents()
				particle_number=recon_event[j].GetPartEventNumber()
				
				tof0_event = spill.GetReconEvents()[j].GetTOFEvent().GetTOFEventSpacePoint().GetTOF0SpacePointArray()
				tof1_event = spill.GetReconEvents()[j].GetTOFEvent().GetTOFEventSpacePoint().GetTOF1SpacePointArray()
				
				# data quality cut: only single space-points in the time-of-flight detectors are allowed. 
				if tof0_event.size() == 1 and tof1_event.size() == 1:
					t0=tof0_event[0].GetTime()
					t1=tof1_event[0].GetTime()
					tdiff01=((t1)-(t0))
					
					# PID cut: only muons are allowed.	
					if tof_min < tdiff01 < tof_max:
						scifi_event = spill.GetReconEvents()[j].GetSciFiEvent()
						if not scifi_event:
							continue
						tracks = scifi_event.scifitracks()
						for k in range(tracks.size()):
							trackpoints = tracks[k].scifitrackpoints()
							for l in range(trackpoints.size()):
								
								# choice of the plane number (0, 1, 2) inside a particular tracker station.
								if trackpoints[l].plane()==plane_number:
									
									# upstream tracker
									if trackpoints[l].tracker()==0:
										# data quality cut: only muons which leave single tracks in upstream tracker are allowed.
										if tracks.size()==1 or tracks.size()==2:
											uspx=trackpoints[l].mom().Px()
											uspy=trackpoints[l].mom().Py()
											uspz=trackpoints[l].mom().Pz()
											usp=np.sqrt(uspx**2+uspy**2+uspz**2)
											if  usp_min < usp < usp_max:
												#p_value = tracks[k].P_value()
												#if p_value > 0.02:
												uspx_cut=trackpoints[l].mom().Px()
												uspy_cut=trackpoints[l].mom().Py()
												uspz_cut=trackpoints[l].mom().Pz()
												usx=trackpoints[l].pos().x()
												usy=trackpoints[l].pos().y()
												usz=trackpoints[l].pos().z()
												station_us=trackpoints[l].station()
												usfile.write('%i %i %i %i %i %g %g %g %g %g %g %g %i %i %i %i %i %i %i %i %i %i %i\n'%(p_0,1,2,0,6-station_us,1,usx,usy,usz,uspx_cut,uspy_cut,uspz_cut,1,1,1,1,0,0,0,0,0,0,0))						
												
									# downstream tracker
									if trackpoints[l].tracker()==1:
										dspx=trackpoints[l].mom().Px()
										dspy=trackpoints[l].mom().Py()
										dspz=trackpoints[l].mom().Pz()
										dsx=trackpoints[l].pos().x()
										dsy=trackpoints[l].pos().y()
										dsz=trackpoints[l].pos().z()
										station_ds=trackpoints[l].station()
										tracks_ds=tracks.size()
										dsfile.write('%i %i %i %i %i %g %g %g %g %g %g %g %i %i %i %i %i %i %i %i %i %i %i\n'%(p_0,1,2,0,station_ds+5,1,dsx,dsy,dsz,dspx,dspy,dspz,1,1,1,1,0,0,0,0,0,0,0))						
										
	return usfile
	return dsfile
	usfile.close()
	dsfile.close()
	
