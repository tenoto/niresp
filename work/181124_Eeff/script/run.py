#!/usr/bin/env python

import os 
import glob

version = 'v181124a'

outdir = 'work/181124_Eeff/effective_area'
cmd = 'rm -rf %s; mkdir -p %s' % (outdir,outdir)
print(cmd);os.system(cmd)

for scaling_factor in [1.00,0.85]:
	for maximum_incident_angle_deg in [3.0,1.4]:
		for foilpara_txt in glob.glob('data/cubesat/foilparameter/v181108/foilpara/*txt'):
			areafile = '%s/%s_i%s_s%s_%s.txt' % (outdir, 
				os.path.splitext(os.path.basename(foilpara_txt))[0],
				str(maximum_incident_angle_deg).replace('.','p'),
				str('%.2f' % scaling_factor).replace('.','p'),
				version)
			cmd  = 'niresp/cli/cubesat.py generateeffectiveareafile %s ' % 'data/cubesat/cubesat.yaml'
			cmd += '--foil_dimension_file %s ' % foilpara_txt
			cmd += '--maximum_incident_angle_deg %.2f ' % maximum_incident_angle_deg
			cmd += '--scaling_factor %.2f ' % scaling_factor			
			cmd += '--outputfile %s' % areafile
			print(cmd);
			try:
				os.system(cmd)
			except:
				print("fail: %s" % areafile)

			if os.path.exists(areafile):
				arffile = areafile.replace('.txt','.arf')
				yamlfile = areafile.replace('.txt','.yaml')				
				cmd  = 'niresp/cli/cubesat.py generatearffile %s ' % yamlfile
				cmd += '%s ' % areafile
				cmd += '%s ' % arffile 
				print(cmd);os.system(cmd)		

rmffile = '%s/gwcubesat_%s.rmf' % (outdir,version)
cmd  = 'niresp/cli/cubesat.py generatermffile %s ' % 'data/cubesat/cubesat.yaml'
cmd += '%s' % rmffile
print(cmd);os.system(cmd)
