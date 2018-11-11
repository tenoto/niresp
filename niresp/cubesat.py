__author__ = 'Teru Enoto'
__version__ = '0.01'
"""
HISTORY

2018-11-10 generated in GitHub repository by Teru Enoto 
"""

import os 
import sys
import yaml
import scipy 
import pandas
import datetime
import numpy as np
import astropy.io.fits as pyfits 

import niresp.xrc 

class CubeSatXRC():
	def __init__(self, yamlfile):
		self.yamlfile = yamlfile 

		print("**** XRC (X-ray Concentrator) created for CubeSat ****")
		print("yamlfile: %s" % self.yamlfile)
		self.readYAMLFile()
		self.setMaterial()		
		self.showProperty()
		self.setFoilDimension()

	def readYAMLFile(self):
		if not os.path.exists(self.yamlfile):
			sys.stderr.write('input parameter_file does not exist.: %s\n' % self.yamlfile)
			quit()		
		f = open(self.yamlfile)
		self.param = yaml.load(f)
		f.close()

		if not (self.param['surface_model'] in ["DW","NC"]):
			print("surface model is not correct.")
			quit()

	def showProperty(self):
		print("------------- Propeties ---------------")
		print("Material : %s" % self.param['material']['symbol'])
		print("   atomicWeight %.3f, Density %.3f (g/cm3)" % (
			self.param['material']['atomicWeight'],
			self.param['material']['density']))
		print("Foil_dimension_file: %s " % self.param['foil_dimension_file'])
		print("Spoke1: number=%d, width=%.3f cm" % (
			self.param['spoke']['n_spoke1'], 
			self.param['spoke']['spoke1_wid']))
		print("Spoke2: number=%d, width=%.3f cm" % (
			self.param['spoke']['n_spoke2'], 
			self.param['spoke']['spoke2_wid']))
		print("Roughness: %.3f (A)" % self.param['roughness_A'])	
		self.use_foilid_list = eval(self.param['use_foilid_list'])
		print("Surface Model: %s" % self.param['surface_model'])

	def setMaterial(self):
		self.material = niresp.xrc.Material(self.param['material']['symbol'], 
			self.param['material']['atomicWeight'],
			self.param['material']['density'])

	def setFoilDimension(self):
		if not os.path.exists(self.yamlfile):
			sys.stderr.write('input parameter_file does not exist.: %s\n' % self.param['foil_dimension_file'])
			quit()		

		self.foils = []
		print("------------------------------------------------------------------------")
		print("Foil Input Dimension")
		print("Foil-ID   : Material, Top radi, Bottom radi, Angle, Angle, GeoArea Use?")
		print("                        (cm)        (cm)     (rad)  (deg)   (cm2)      ")
		print("------------------------------------------------------------------------")
		total_geometric_area     = 0.0
		total_geometric_area_obs1 = 0.0
		total_geometric_area_obs2 = 0.0		
		for line in open(self.param['foil_dimension_file']):
			cols = line.split() 
			if cols[0] == '#':
				continue 

			foil_id = len(self.foils)+1

			foil_parab_top_radius = float(cols[0]) # mm
			foil_parab_bot_radius = float(cols[1]) # mm 
			foil_parab_angle_rad  = float(cols[2])
			foil_parab_angle_deg  = scipy.degrees(foil_parab_angle_rad)

			if foil_parab_angle_deg > self.param['maximum_incident_angle_deg']:
				break 

			foil_parab_bot_radius_cm = 0.1 * foil_parab_bot_radius
			foil_parab_top_radius_cm = 0.1 * foil_parab_top_radius

			foil_cone_angle_rad = foil_parab_angle_rad # assumed for CubeSat 

			self.foils.append(niresp.xrc.Foil(foil_id, 
				self.material, 
				foil_parab_top_radius_cm, foil_parab_bot_radius_cm, 
				foil_cone_angle_rad))

			self.foils[-1].roughness_A = self.param['roughness_A']
			self.foils[-1].n_spoke1    = self.param['spoke']['n_spoke1']
			self.foils[-1].n_spoke2    = self.param['spoke']['n_spoke2']
			self.foils[-1].spoke1_wid  = self.param['spoke']['spoke1_wid']
			self.foils[-1].spoke2_wid  = self.param['spoke']['spoke2_wid']
			self.foils[-1].setUseFlag(self.use_foilid_list)
			self.foils[-1].surface_model = self.param['surface_model']
			self.foils[-1].showProperty()

			foil_geometric_area = (foil_parab_top_radius_cm**2-foil_parab_bot_radius_cm**2) * np.pi

			total_geometric_area      += foil_geometric_area
			total_geometric_area_obs1 += foil_geometric_area * self.foils[-1].getSpokeObsculation()
			subtraction  = float(self.param['spoke']['n_spoke1']) * float(self.param['spoke']['spoke1_wid'])*10.0 * (self.foils[-1].top_radius_cm-self.foils[-1].bottom_radius_cm)*10.0
			subtraction += float(self.param['spoke']['n_spoke2']) * float(self.param['spoke']['spoke2_wid'])*10.0 * (self.foils[-1].top_radius_cm-self.foils[-1].bottom_radius_cm)*10.0
			value = foil_geometric_area - subtraction
			total_geometric_area_obs2 += value

		print("Total geometric area (cm2):",        total_geometric_area)
		print("after subtracting spoke obsculation: ", total_geometric_area_obs2)
		obst_1 = float(self.param['spoke']['n_spoke1']) * float(self.param['spoke']['spoke1_wid']) * (self.foils[-1].top_radius_cm-self.foils[0].bottom_radius_cm)
		obst_2 = float(self.param['spoke']['n_spoke2']) * float(self.param['spoke']['spoke2_wid']) * (self.foils[-1].top_radius_cm-self.foils[0].bottom_radius_cm)
		print("------------------------------------------------------------------------")

	def getFoil(self, foilid):
		return self.foils[foilid-1]

	def set_foil_dimension_file(self,foil_dimension_file):
		self.param['foil_dimension_file'] = foil_dimension_file

	def set_maximum_incident_angle_deg(self,maximum_incident_angle_deg):
		self.param['maximum_incident_angle_deg'] = maximum_incident_angle_deg

	def getEffectiveArea(self, energy_keV):
		Aeff_tot = 0.0
		for foil in self.foils:
			Aeff_tot += foil.getEffectiveArea(energy_keV)
		return scipy.real(Aeff_tot)

	def generateEffectiveAreaFile(self,outfile='effectivearea.txt'):
		print("generateEffectiveAreaFile")
		energy_array = np.arange(
			self.param['effectivearea_energy_start_keV'],
			self.param['effectivearea_energy_stop_keV'],
			self.param['effectivearea_energy_step_keV'])
		f = open(outfile,'w')		
		dump = '# foil numbers : %d\n' % len(self.foils)
		f.write(dump)
		for energy_keV in energy_array:
			dump = '%.4f %.4f\n' % (energy_keV,self.getEffectiveArea(energy_keV))
			f.write(dump)
		f.close()

#	def generateArffile(self,outfile='gwcubesat.arf'):
#		print("generateArffile")


class CubeSatXspecResponse():
	def __init__(self,setup_yamlfile_path):
		self.setup_yamlfile_path = setup_yamlfile_path

		if not os.path.exists(self.setup_yamlfile_path):
			raise FileNotFoundError("{} not found".format(self.setup_yamlfile_path))
		try:
			self.param = yaml.load(open(self.setup_yamlfile_path))
		except OSError as e:
			raise		
		print("XspecResponse: parameter file {} is loaded.".format(self.setup_yamlfile_path))

		self.NENERGY = int((self.param['ENERGY_MAX']-self.param['ENERGY_MIN'])/self.param['ENERGY_STEP']) # number of column in EBOUNDS 
		self.NPI     = int((self.param['PI_MAX']-self.param['PI_MIN'])/self.param['PI_STEP'])
		self.DETCHANS = self.NPI+1
		print('---- NENERGY={}, NPI={}, DETCHANS={}'.format(self.NENERGY,self.NPI,self.DETCHANS))

		# 1st extensiton of arf file and 2nd extension of rmffile 
		self.energy_lo_array = np.array([self.param['ENERGY_MIN']+self.param['ENERGY_STEP']*i for i in range(0,self.NENERGY+1)])
		self.energy_hi_array = np.array([self.param['ENERGY_MIN']+self.param['ENERGY_STEP']*(i+1) for i in range(0,self.NENERGY+1)])	
		self.energy_cen_array = 0.5*(self.energy_lo_array+self.energy_hi_array)	

	def generate_rmffile(self,outrmffile='gwcubesat.rmf'):
		"""
		this method generates XPSEC redistribution matrix file (rmf file) for NICER.
		https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html#tth_sEc3.2
		"""		
		print("--------------------------------------------")
		print("NICER rmf generator (ver.{})".format(__version__))
		print(" output file : {} ".format(outrmffile))
		print("--------------------------------------------")

		outdir = os.path.dirname(outrmffile)
		if outdir == '':
			outdir = './'
		elif not os.path.exists(outdir):
			os.makedirs(outdir)
		if os.path.exists(outrmffile):
			sys.stdout.write("{} already exists.\n".format(outrmffile))
			quit()

		# --- prepare EBOUNDS (1st Extension) --- 
		channel_array = np.array([i for i in range(0,self.NPI+1)])
		e_min_array = np.array([self.param['PI_MIN']+self.param['PI_STEP']*i for i in range(0,self.NPI+1)])
		e_max_array = np.array([self.param['PI_MIN']+self.param['PI_STEP']*(i+1) for i in range(0,self.NPI+1)])	
		e_cen_array = 0.5*(e_min_array+e_max_array)	
		e_min_column = pyfits.Column(name='E_MIN',format="E",array=e_min_array,unit="keV")
		e_max_column = pyfits.Column(name='E_MAX',format="E",array=e_max_array,unit="keV")	
		channel_column = pyfits.Column(name='CHANNEL',format='J',array=channel_array,unit="chan")
		hdu_ebounds = pyfits.ColDefs([channel_column,e_min_column,e_max_column])

		# --- read ascii rmf file 
		rmf_table = self.param['RMF_TABLE']
		if not os.path.exists(rmf_table):
			print("Error: rmf ascii file {} does not exist.".format(rmf_table))
			quit()
		print("-- start reading ascii rmf file: {}".format(rmf_table))

		numoflines = len(open(rmf_table).readlines())
		rmf_2darray = np.zeros((self.NENERGY+1,self.NPI+1))
		i = 0
		for line in open(rmf_table):	
			cols = line.split() 
			if len(cols) != 5 or cols[0] == 'i':
				continue 
			ien = int(cols[0])   # the energy grid index (also starting at 0)
			iph = int(cols[1])   # the pulse-height index  (starting at 0)
			en  = float(cols[2]) # the central bin energy in eV
			ph  = float(cols[3]) # the pulse-height energy in eV
			val = float(cols[4]) # the normalized redistribution element
			#acut = float(cols[5]) # normalized entry when the cutoff is included 
			#rmf_2darray[ien][iph] = val
			rmfval  = val 
			rmf_2darray[ien][iph] = rmfval 
			if i % 200000 == 0:
				sys.stdout.write('...%.1f%% processed\n' % (100.0*float(i)/float(numoflines)))		
			i += 1
		print("-- finished reading ascii rmf file.")	

		matrix_form = "%sE" % (self.NPI+1)
		matrix_column = pyfits.Column(name='MATRIX',format=matrix_form, array=rmf_2darray)
		#matrix_column = pyfits.Column(name='MATRIX',format=matrix_form, array=rmf_2darray_norm)
		energy_lo_column = pyfits.Column(name='ENERG_LO',format="E",array=self.energy_lo_array,unit="keV")
		energy_hi_column = pyfits.Column(name='ENERG_HI',format="E",array=self.energy_hi_array,unit="keV")	
		n_grp_column = pyfits.Column(name='N_GRP',format="J",array=np.array([self.param['N_GRP'] for i in range(0,self.NENERGY+1)]),unit="")
		f_chan_column = pyfits.Column(name='F_CHAN',format="1J",array=np.array([self.param['F_CHAN'] for i in range(0,self.NENERGY+1)]),unit="")	
		n_chan_column = pyfits.Column(name='N_CHAN',format="1J",array=np.array([self.NPI+1 for i in range(0,self.NENERGY+1)]),unit="")			
		hdu_matrix = pyfits.ColDefs([energy_lo_column,energy_hi_column,n_grp_column,f_chan_column,n_chan_column,matrix_column])

		# --- Filling the extensions to a fits file ---
		prhdu = pyfits.PrimaryHDU()
		tbhdu_ebounds = pyfits.BinTableHDU.from_columns(hdu_ebounds)
		tbhdu_ebounds.name = 'EBOUNDS'
		tbhdu_matrix = pyfits.BinTableHDU.from_columns(hdu_matrix)
		tbhdu_matrix.name = 'SPECRESP MATRIX'	
		hdu = pyfits.HDUList([prhdu,tbhdu_ebounds,tbhdu_matrix])

		# --- Write to output fitsfile --- 
		hdu.writeto(outrmffile)
		print("-- write to fits file: {}".format(outrmffile))

		f = open('%s/tmp_header_cmn.txt' % outdir,'w')
		dump  = "COMMENT gennerated by script %s (version %s)\n" % (sys.argv[0],__version__)
		dump += "COMMENT at %s\n" % datetime.datetime.now()
		dump += "COMMENT parameter file : %s\n" % self.setup_yamlfile_path
		#for key in param:
		#	dump += "COMMENT %s: %s\n" % (key, param[key])
		for line in open(self.setup_yamlfile_path):
			dump += "COMMENT %s\n" % line 
		f.write(dump)
		f.close()


		f = open('%s/tmp_header_ext1.txt' % outdir,'w')
		dump = """HDUCLASS= 'OGIP    '           / format conforms to OGIP standard
HDUCLAS1= 'RESPONSE'           / dataset relates to spectral response
HDUCLAS2= 'EBOUNDS '           / nominal energies of PHA chan boundaries
HDUVERS = '1.3.0   '           / Version of format (OGIP memo CAL/GEN/92-002a)
HDUDOC  = 'OGIP memos CAL/GEN/92-002 & 92-002a' / Documents describing the forma
HDUVERS1= '1.0.0   '           / Obsolete - included for backwards compatibility
HDUVERS2= '1.3.0   '           / Obsolete - included for backwards compatibility
TELESCOP= 'NICER   '           / mission/satellite name
INSTRUME= 'XTI     '           / instrument/detector name
FILTER  = 'XTI     '           / filter in use
DETCHANS=                 %d / total number of detector channels
CHANTYPE= 'PI      '           / WARNING This is NOT an OGIP-approved value
RMFVERSN= '1992a   '           / Obsolete - included for backwards compatibility
TLMIN1  =                   %d / Minimum value legally allowed in column 1
TLMAX1  =                 %d / Maximum value legally allowed in column 1
HISTORY EBOUNDS extension written by wtebd3 1.1.1	
"""	% (self.DETCHANS,channel_array[0],channel_array[-1])
		f.write(dump)
		f.close()

		f = open('%s/tmp_header_ext2.txt' % outdir,'w')
		dump = """EXTNAME = 'SPECRESP MATRIX'    / name of this binary table extension
HDUCLASS= 'OGIP    '           / format conforms to OGIP standard
HDUCLAS1= 'RESPONSE'           / dataset relates to spectral response
HDUVERS1= '1.0.0   '           / version of family of formats
HDUCLAS2= 'RSP_MATRIX'         / datasets is a spectral response matrix
HDUVERS2= '1.2.0   '           / Version of format (OGIP memo CAL/GEN/92-002a)
HDUCLAS3= 'FULL    '           / includes all efficiencies
TELESCOP= 'NICER   '           / mission/satellite
INSTRUME= 'XTI     '           / mission/satellite
DETNAM  = '        '           / detector in use
FILTER  = 'XTI     '           / filter in use
DETCHANS=                 %d / total number of detector channels
CHANTYPE= 'PI      '           / detector channel type (PHA or PI)
HIERARCH LO_THRESH =  %.2e / lower threshold for stored matrix
LO_THRES=             %.2e / lower threshold for stored matrix
RMFVERSN= '1992a   '           / deprecated classification of format
TLMIN4  =                    %d / Start channel used
TLMAX4  =                 %d / Last channel used	
"""	% (self.DETCHANS, self.param['LO_THRES'], self.param['LO_THRES'], channel_array[0],channel_array[-1])
		f.write(dump)
		f.close()	

		print('-- editing the header keywords...')

		cmd  = 'fthedit %s+1 @%s/tmp_header_ext1.txt;\n' % (outrmffile,outdir)
		cmd += 'fthedit %s+2 @%s/tmp_header_ext2.txt;\n' % (outrmffile,outdir)	
		for i in range(3):
			cmd += 'fthedit %s+%d @%s/tmp_header_cmn.txt;\n' % (outrmffile,i,outdir)	
			cmd += 'fparkey %s %s+%d NIRESPVER add=yes;\n' % (__version__,outrmffile,i)
		print(cmd); os.system(cmd)
		os.system('rm -f %s/tmp_header_*.txt' % outdir)
		print('DONE!\n')

	def generate_arffile(self,inareafile,outarffile):
		"""
		this method generates XPSEC redistribution matrix file (rmf file) for NICER.
		https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html#tth_sEc3.2
		"""		
		print("--------------------------------------------")
		print("NICER arf generator (ver.{})".format(__version__))
		print(" output file : {} ".format(outarffile))		
		print("--------------------------------------------")

		outdir = os.path.dirname(outarffile)
		if outdir == '':
			outdir = './'
		elif not os.path.exists(outdir):
			os.makedirs(outdir)
		if os.path.exists(outarffile):
			sys.stdout.write("{} already exists.\n".format(outarffile))
			quit()

		energy_lo_column = pyfits.Column(name='ENERG_LO',format="E",array=self.energy_lo_array,unit="keV")
		energy_hi_column = pyfits.Column(name='ENERG_HI',format="E",array=self.energy_hi_array,unit="keV")	

		"""
		xrc_snid_table = self.param['XRC_SNID_TABLE']
		xrc_list = []
		for line in open(xrc_snid_table):
			cols = line.split()
			if cols[0] == '#':
				continue 
			sn = cols[3]
			science_detector_id = cols[4]
			mpu_num = science_detector_id[0]
			fpm_num = science_detector_id[1]
			sys.stdout.write('=======================================\n')
			sys.stdout.write('... seting SN=%s SciID=%s MPU=%s FPM=%s\n' % 
				(sn, science_detector_id, mpu_num, fpm_num))

			if int(sn) in self.param['EXCLUDE_SN']:
				sys.stdout.write('...skipping SN=%s\n' % sn)				
				continue 

			xrc_list.append(NicerSingleModule(sn))
			xrc_list[-1].setMPUnumber(mpu_num)
			xrc_list[-1].setFPMnumber(fpm_num)		
			xrc_list[-1].readEffectiveArea(self.param['XRC_AREA_DIR'])
			xrc_list[-1].readApertureFraction(self.param['XRC_APERTURE_TABLE'])
			xrc_list[-1].area_cm2_array = np.array([np.interp(energy_cen,xrc_list[-1].energy_keV_list,xrc_list[-1].area_cm2_aperture_list) 
					for energy_cen in self.energy_cen_array])	

		allxrc_area_cm2_array = np.array([0.0 for energy_cen in self.energy_cen_array])
		for xrc in xrc_list:
			sys.stdout.write('...adding effective area of SN=%s\n' % xrc.SN)			
			allxrc_area_cm2_array += xrc.area_cm2_array
		"""

		df_area = pandas.read_table(inareafile,
			header=None,names=('energy_keV','area_cm2'),
			index_col=None,sep=' ',skiprows=1)

		allxrc_area_cm2_array = np.array(
			[self.param['XRC_APERTURE_FRACTION']/100.0 * np.interp(energy_cen,df_area['energy_keV'],df_area['area_cm2']) 
			for energy_cen in self.energy_cen_array])

		# --- set QE ---
		sys.stdout.write('reading QE table.\n')
		if not os.path.exists(self.param['QE_TABLE']):
			sys.stderr.write("QE table %s does not exist.\n" % self.param['QE_TABLE'])
			quit()
		qe_energy_keV_list = []
		qe_efficency_list = []
		for line in open(self.param['QE_TABLE']):
			cols = line.split()
			if cols[0] == 'energy_eV':
				continue 
			energy_eV = float(cols[0])
			qe_without_trigg_efficiency = float(cols[1])
			qe_energy_keV_list.append(energy_eV*1e-3)
			qe_efficency_list.append(qe_without_trigg_efficiency)	
		qe_efficency_array = np.array([np.interp(energy_cen,qe_energy_keV_list,qe_efficency_list) for energy_cen in self.energy_cen_array])

		# --- set window transmission table --- 
		sys.stdout.write('reading window transmission table.\n')
		if not os.path.exists(self.param['WINDOW_TABLE']):
			sys.stderr.write("Window transmission table %s does not exist.\n" % self.param['WINDOW_TABLE'])
			quit()
		window_energy_keV_list = []
		window_efficiency_list = []
		for line in open(self.param['WINDOW_TABLE']):
			cols = line.split()
			if cols[0] == 'Energy_eV':
				continue 
			energy_eV = float(cols[0])
			window_transmission = float(cols[1])
			window_energy_keV_list.append(energy_eV*1e-3)
			window_efficiency_list.append(window_transmission)	
		window_efficiency_array = np.array([np.interp(energy_cen,window_energy_keV_list,window_efficiency_list) for energy_cen in self.energy_cen_array])

		# --- Thermal shield ---
		sys.stdout.write('reading thermal shield table file.\n')
		if not os.path.exists(self.param['THERMALSD_TABLE']):
			sys.stderr.write("Thermal shield table %s does not exist.\n" % self.param['THERMALSD_TABLE'])
			quit()
		thermalsd_energy_keV_list = []
		thermalsd_transmission_list = []
		for line in open(self.param['THERMALSD_TABLE']):
			cols = line.split()
			if cols[0] == 'SX700-Energy':
				continue 
			energy_eV = float(cols[0])
			thermalsd_transmission = float(cols[1])
			thermalsd_energy_keV_list.append(energy_eV*1e-3)
			thermalsd_transmission_list.append(thermalsd_transmission)	
		thermalsd_transmission_array = np.array([np.interp(energy_cen,thermalsd_energy_keV_list,thermalsd_transmission_list)  for energy_cen in self.energy_cen_array])

		total_area_cm2_array = allxrc_area_cm2_array * qe_efficency_array * window_efficiency_array * thermalsd_transmission_array

		ecen_column = pyfits.Column(name='ENERGY',format='E',array=self.energy_cen_array,unit="keV")
		allxrc_area_cm2_column = pyfits.Column(name='XRCAREA',format='E',array=allxrc_area_cm2_array,unit="cm2")
		qe_efficiency_column = pyfits.Column(name='QE',format="E",array=qe_efficency_array,unit="")
		window_efficiency_column = pyfits.Column(name='WINDOW',format="E",array=window_efficiency_array,unit="")	
		thermalsd_transmission_column = pyfits.Column(name='THERMALSD',format="E",array=thermalsd_transmission_array,unit="")

		specresp_column = pyfits.Column(name='SPECRESP',format='E',array=total_area_cm2_array,unit="cm**2")
		hdu_specresp = pyfits.ColDefs([
			energy_lo_column,
			energy_hi_column,
			specresp_column,
			ecen_column,
			allxrc_area_cm2_column,
			qe_efficiency_column,
			window_efficiency_column,
			thermalsd_transmission_column])

		# --- Filling the extensions to a fits file ---
		prhdu = pyfits.PrimaryHDU()
		tbhdu_specresp = pyfits.BinTableHDU.from_columns(hdu_specresp)
		tbhdu_specresp.name = 'SPECRESP'
		hdu = pyfits.HDUList([prhdu,tbhdu_specresp])

		# --- Write to output fitsfile --- 
		hdu.writeto(outarffile)

		f = open('%s/tmp_header_cmn.txt' % outdir,'w')
		dump  = "COMMENT gennerated by script xspecresp.py (version %s)\n" % (__version__)
		dump += "COMMENT at %s\n" % datetime.datetime.now()
		dump += "COMMENT parameter file : %s\n" % self.setup_yamlfile_path
		for line in open(self.setup_yamlfile_path):
			dump += "COMMENT %s\n" % line 
		f.write(dump)
		f.close()

		f = open('%s/tmp_header_specresp.txt' % outdir,'w')
		dump = """HDUCLASS= 'OGIP    '           / format conforms to OGIP standard
HDUCLAS1= 'RESPONSE'           / dataset relates to spectral response
HDUCLAS2= 'SPECRESP'           / dataset contains spectral response
HDUVERS = '1.1.0   '           / Version of format (OGIP memo CAL/GEN/92-002a)
HDUDOC  = 'OGIP memos CAL/GEN/92-002 & 92-002a' / Documents describing the forma
HDUVERS1= '1.0.0   '           / Obsolete - included for backwards compatibility
HDUVERS2= '1.1.0   '           / Obsolete - included for backwards compatibility
TELESCOP= 'NICER   '           / mission/satellite name
INSTRUME= 'XTI     '           / instrument/detector name
DETNAM  = '        '           / specific detector name in use
FILTER  = 'XTI     '           / filter in use
ARFVERSN= '1992a   '           / Obsolete - included for backwards compatibility
PHAFILE = 'UNKNOWN '           / PHA file for which this ARF created
HISTORY This file is made with the ftool addarf 1.2.7
HISTORY  FITS ARF extension written by WTARF1 1.1.0
"""	
		f.write(dump)
		f.close()

		cmd = 'fthedit %s+1 @%s/tmp_header_specresp.txt;' % (outarffile,outdir)
		for i in range(2):
			cmd += 'fthedit %s+%d @%s/tmp_header_cmn.txt;' % (outarffile,i,outdir)	
		for i in range(2):
			cmd += 'fparkey %s %s+%d NIRESPVER add=yes;' % (__version__,outarffile,i)
		print(cmd);os.system(cmd)

		os.system('rm -f %s/tmp_header_*.txt' % outdir)

class NicerSingleModule:
	def __init__(self,SN):
		self.SN = SN 
		sys.stdout.write('---- reading SN%s...\n' % self.SN)

	def setMPUnumber(self, mpu_num):
		self.mpu_num = mpu_num
		sys.stdout.write('------ set mpu_num = %s\n' % mpu_num)

	def setFPMnumber(self, fpm_num):
		self.fpm_num = fpm_num
		sys.stdout.write('------ set fpm_num = %s\n' % fpm_num)

	def readEffectiveArea(self,directory):
		eff_filename = '%s/xrc%s_area.txt' % (directory,self.SN)
		
		self.energy_keV_list = []
		self.area_cm2_list   = []
		for line in open(eff_filename):
			cols = line.split() 
			energy_keV = float(cols[0])
			area_cm2   = float(cols[1])
			self.energy_keV_list.append(energy_keV)
			self.area_cm2_list.append(area_cm2)

	def readApertureFraction(self,aperturefile):
		for line in open(aperturefile):
			cols = line.split()
			if cols[0] == '#S/N':
				continue
			if int(cols[0]) == int(self.SN):
				self.aperture_fraction = 0.01 * float(cols[2])
				break
		self.area_cm2_aperture_list = [self.aperture_fraction * i for i in self.area_cm2_list]

	"""
	2018-02-10: Trigger efficiency is moved to the rmf file, already included in the input table generated by Jack.

	def readTriggerEfficiency(self,trigefffile):
		self.error_func_center_pi = None
		self.error_func_center_pi_error = None
		self.error_func_sigma_pi = None
		self.error_func_sigma_pi_error = None 
		for line in open(trigefffile):
			cols = line.split()
			if len(cols) == 0:
				continue 
			if cols[0] in ['Error','MPU']:
				continue
			trig_mpu_num = cols[0]
			trig_fpm_num = cols[1]			
			if (int(trig_mpu_num)==int(self.mpu_num)) and (int(trig_fpm_num)==int(self.fpm_num)):
				self.error_func_center_pi = float(cols[2])
				self.error_func_center_pi_error = float(cols[3])
				self.error_func_sigma_pi = float(cols[4])
				self.error_func_sigma_pi_error = float(cols[5])	
				break 

		sys.stdout.write('trig_mpu_num=%s trig_fpu_num=%s error function: center_pi=%s (error %s) sigma_pi=%s (error %s) \n' % (
			trig_mpu_num,trig_fpm_num,
			self.error_func_center_pi,self.error_func_center_pi_error,
			self.error_func_sigma_pi,self.error_func_sigma_pi_error))		
		self.thresh_ecen_keV  = niPI2keV(self.error_func_center_pi) # centroid of efficiency function (keV)
		self.thresh_sigma_keV = niPI2keV(self.error_func_sigma_pi) # sigma width of efficiency function (keV)
		sys.stdout.write('-- error function: center_keV=%s sigma_keV=%s\n' % (
			self.thresh_ecen_keV,self.thresh_sigma_keV))

		self.trig_efficiency_list = []
		self.area_cm2_aperture_trig_list = []
		for i in range(len(self.energy_keV_list)):
			e_cen = self.energy_keV_list[i]
			trig_eff_factor = 0.5 * (1.0+special.erf((e_cen-self.thresh_ecen_keV)/(self.thresh_sigma_keV*numpy.sqrt(2.0))))
			self.area_cm2_aperture_trig_list.append(trig_eff_factor * self.area_cm2_aperture_list[i])
			self.trig_efficiency_list.append(trig_eff_factor)
	"""




