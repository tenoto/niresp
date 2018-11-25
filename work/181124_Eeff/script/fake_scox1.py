#!/usr/bin/env python

import os 
import glob
import astropy.io.fits as fits

def run_fake(arffile,exposure_ks):
	outdir = '%s/fake_scox1' % os.path.dirname(arffile)
	basename = os.path.basename(arffile)
	outpha = '%s/%s_%sks.pha' % (outdir,os.path.splitext(basename)[0],str(exposure_ks).replace('.','p'))
	outxcm = '%s/%s_%sks.xcm' % (outdir,os.path.splitext(basename)[0],str(exposure_ks).replace('.','p'))	
	outlog = '%s/%s_%sks.log' % (outdir,os.path.splitext(basename)[0],str(exposure_ks).replace('.','p'))
	outps= '%s/%s_%sks.ps' % (outdir,os.path.splitext(basename)[0],str(exposure_ks).replace('.','p'))			
	outpdf= '%s/%s_%sks.pdf' % (outdir,os.path.splitext(basename)[0],str(exposure_ks).replace('.','p'))				
	cmd = 'rm -f %s %s %s %s %s' % (outpha,outxcm,outlog,outps,outpdf)
	print(cmd);os.system(cmd)
	cmd = 'mkdir -p %s' % outdir
	print(cmd);os.system(cmd)

	focal_length_cm = (float(arffile.split('_fl')[-1].split('cm')[0].replace('p','.')))

	hdu = fits.open(arffile)
	maximum_diameter_cm = hdu['SPECRESP'].header['MAXDIAMT']		
	maximum_incident_angle_deg = hdu['SPECRESP'].header['MAXANGL']	
	num_of_foil = hdu['SPECRESP'].header['NUMFOILS']
	scaleing_factor = hdu['SPECRESP'].header['SCALINGF']

	cmd  = "xspec<<EOF\n"
	cmd += 'data 1 data/cubesat/backgrnd/single_xrc_backgrnd.pha\n'
	cmd += 'back 1 data/cubesat/backgrnd/single_xrc_backgrnd.pha\n'	
	cmd += "@work/181124_Eeff/script/ScoX1_Revnivtsev2014_Fig7_approximated.xcm\n"
	cmd += "log %s\n" % outlog
	cmd += "fakeit none\n"
	cmd += "work/181124_Eeff/effective_area/gwcubesat_v181124a.rmf\n"
	cmd += "%s\n" % arffile
	cmd += "y\n"
	cmd += "\n"
	cmd += "%s\n" % outpha
	cmd += "%d\n" % (exposure_ks*1000.0)
	cmd += "setplot energy\n"
	cmd += "setplot rebin 3 100\n"
	cmd += "ignore **-0.2 3.0-**\n"
	cmd += "show rate\n"
	cmd += "cpd /xw\n"
	cmd += "plot ld\n"
	cmd += "ipl ld\n"
	cmd += "lwid 5\n"
	cmd += "lwid 5 on 1..100\n"
	cmd += "col 2 on 2\n"
	cmd += "time off\n"
	cmd += "r x 0.2 3.0\n"
	cmd += "r y 0.3 300.0\n"
	cmd += "la t %s\n" % os.path.basename(outpha)
	cmd += "hard %s/cps\n" % outps
	cmd += "exit\n"
	cmd += "save all %s\n" % outxcm
	cmd += "log none\n"
	cmd += "exit\n"
	cmd += "EOF\n"
	print(cmd);os.system(cmd)	

	cmd = 'ps2pdf %s' % outps
	print(cmd);os.system(cmd)	

	cmd = 'mv %s work/181124_Eeff/effective_area/fake_scox1' % os.path.basename(outpdf)
	print(cmd);os.system(cmd)	

	cmd = 'rm -f %s' % (outps)
	print(cmd);os.system(cmd)

	for line in open(outlog):
		cols = line.split()
		if cols[0] == '#Net':
			count_rate = float(cols[6])
			break
	return focal_length_cm, count_rate, num_of_foil, maximum_incident_angle_deg, maximum_diameter_cm, scaleing_factor

exposure_ks = 10.0
f = open('work/181124_Eeff/effective_area/fl_vs_rate_scox1.txt','w')
f.write('focal_length_cm, count_rate, num_of_foil, maximum_incident_angle_deg, maximum_diameter_cm, scaleing_factor \n')
for arffile in glob.glob('work/181124_Eeff/effective_area/*.arf'):
	print(arffile)
	f.write('%.2f %.2f %d %.2f %.2f %.3f \n' % run_fake(arffile,exposure_ks))
f.close()

exposure_ks = 0.06
run_fake('work/181124_Eeff/effective_area/foilpara_fl24p92cm_i3p0_s1p00_v181124a.arf',exposure_ks)
run_fake('work/181124_Eeff/effective_area/foilpara_fl50cm_i3p0_s1p00_v181124a.arf',exposure_ks)
exit()

exposure_ks = 1.00
run_fake('work/181124_Eeff/effective_area/foilpara_fl24p92cm_i3p0_s1p00_v181124a.arf',exposure_ks)
run_fake('work/181124_Eeff/effective_area/foilpara_fl50cm_i3p0_s1p00_v181124a.arf',exposure_ks)


