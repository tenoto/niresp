# niresp
NICER Response Generator

## 1. Scripts and Data 

Python script for our XRC semi-analytical model, and response generators of rmf and arf files:

	https://github.com/tenoto/niresp

Data files to be used for the above calculation (Since these are NICER internal information and not open to public)

	[xxx] please ask Teru
	# Please download the “data” directory to under “niresp” directory

## 2. How to setup

I assumed …
	- “GitHub” and “pipenv” are available under your environment.

You will be ready to use script after using following command lines≥

	git clone https://github.com/tenoto/niresp.git
	cd niresp 
	cp -r data niresp 
	pipenv install

This command install required python libraries to your environment with making virtual python environment. 
	
	source setenv/setenv.bashrc 

This sets up environmental parameters.

## 3. Script definitions 

	niresp/xrc.py : This library calculates semi-analytical effective area, and makes output files. 

	niresp/xspecresp.py : This convert input text format files to “rmf” and “arf” file for XSPEC analyses of NICER

	niresp/cli/*.py : These files support you to use the above library in command lines.

## 4. How to use 

The following command line gives you which sub command lines exist. 

	niresp/cli/xspecresp.py —help

You can also check sub-command line as :

	niresp/cli/xspecresp.py generate-rmffile —help

If you want to generate “rmf” and “arf” files, following shell script is a test file for this purpose. 

	tests/test_xspecresp.sh

(Yaml parameters are just sample) 

## 5. Note

- The present version is not the final version.  Purpose of this distribution is just sharing the code and to double-check whether we do not have miscommunication of definitions of data and script. This version only shows a framework and code for my semi-analytical calculation and response generations. There are still several points which I would like to improve…. 

- Parameters are tentative. I need to more carefully check parameters (surface roughness, density) used in the script because there have been continuing several discussions and trial to find out the best parameter set. 

- original direcotry : /Users/enoto/Dropbox/enoto/library/xrcAeffCalc/xtisim_1.05
