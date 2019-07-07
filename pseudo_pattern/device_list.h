#ifndef DEVICE_LIST_H
#define DEVICE_LIST_H


TString detector_channel[]={"sam1","sam3","sam5",  	 // index; 0-2
			    "sam2","sam4","sam6","sam8", // index:  3-6
			    "usl","usr","dsl","dsr",	 // index : 7 - 10
			    "bcm_an_ds3"};

TString monitor_channel[]={"bpm4aX","bpm4aY",
			   "bpm4eX","bpm4eY",
			   "bpm12X"};

TString combiner_channel[]={"asym_us_avg","asym_us_dd",
			    "asym_ds_avg","asym_ds_dd",
			    "asym_left_avg","asym_left_dd",
			    "asym_right_avg","asym_right_dd",
			    "asym_sam_15_avg","asym_sam_15_dd"};

Double_t weight[10][2]={ {0.5,0.5}, {0.5,-0.5},
			 {0.5,0.5}, {0.5,-0.5},
			 {0.5,0.5}, {0.5,-0.5},
			 {0.5,0.5}, {0.5,-0.5},
			 {0.5,0.5}, {0.5,-0.5}};
// Uhhh......
Int_t det_pair[10][2] = { {7,8},{7,8},
			  {9,10},{9,10},
			  {7,9},{7,9},
			  {8,10},{8,10},
			  {0,2},{0,2}};
		       


#endif
