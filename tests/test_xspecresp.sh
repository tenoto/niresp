#!/bin/sh -f 

niresp/cli/xspecresp.py generate-rmffile \
	data/xspecresp/xspecresp_setup_param_v181105.yaml \
	out/nicer_v1.04.rmf

niresp/cli/xspecresp.py generate-arffile \
	data/xspecresp/xspecresp_setup_param_v181105.yaml \
	out/ni_xrcall_onaxis_v1.04.arf 
