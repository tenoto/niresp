__author__ = 'Teru Enoto'
__version__ = '1.06'
"""
HISTORY

2018-11-05 generated in GitHub repository by Teru Enoto from xtisim 1.05
"""

#class XrayConcentrator():

# --------------------------------------------
# Libraries
# --------------------------------------------
import os 
import sys 
import yaml
import scipy
from scipy  import interpolate
from scipy.constants import *
Re_cm = physical_constants['classical electron radius'][0] * 1e+2 # classical electron radius 2.818e-13  (cm)
A_cm  = scipy.constants.angstrom * 1e+2 # angstrom in cm unit
ReflectionIndexVacuum = 1.0

# --------------------------------------------
# Common Function
# --------------------------------------------
def convertEnergykeV2LambdaCM(energy_keV):
	lambda_cm = 1.240e-7 / energy_keV
	return lambda_cm

def convertEnergykeV2LambdaA(energy_keV):
	lambda_A = 12.39 / energy_keV # Angstrom
	return lambda_A

# --------------------------------------------
# Main Body 
# --------------------------------------------
class XrayConcentrator:
	def __init__(self, yamlfile):
		self.yamlfile = yamlfile 

		print "**** XRC (X-ray Concentrator) created ****"
		print "yamlfile: %s" % self.yamlfile
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
			print "surface model is not correct."
			quit()

	def showProperty(self):
		print "------------- Propeties ---------------"
		print "Material : %s" % self.param['material']['symbol']
		print "   atomicWeight %.3f, Density %.3f (g/cm3)" % (
			self.param['material']['atomicWeight'],
			self.param['material']['density'])
		print "Foil_dimension_file: %s " % self.param['foil_dimension_file']
		print "Spoke1: number=%d, width=%.3f cm" % (
			self.param['spoke']['n_spoke1'], 
			self.param['spoke']['spoke1_wid'])
		print "Spoke2: number=%d, width=%.3f cm" % (
			self.param['spoke']['n_spoke2'], 
			self.param['spoke']['spoke2_wid'])		
		print "Roughness: %.3f (A)" % self.param['roughness_A']		
		self.use_foilid_list = eval(self.param['use_foilid_list'])
		print "Surface Model: %s" % self.param['surface_model']

	def setMaterial(self):
		self.material = Material(self.param['material']['symbol'], 
			self.param['material']['atomicWeight'],
			self.param['material']['density'])

	def setFoilDimension(self):
		if not os.path.exists(self.yamlfile):
			sys.stderr.write('input parameter_file does not exist.: %s\n' % self.param['foil_dimension_file'])
			quit()		

		self.foils = []
		print "------------------------------------------------------------------------"
		print "Foil Input Dimension"
		print "Foil-ID   : Material, Top radi, Bottom radi, Angle, Angle, GeoArea Use?"
		print "                        (cm)        (cm)     (rad)  (deg)   (cm2)      "
		print "------------------------------------------------------------------------"
		total_geometric_area     = 0.0
		total_geometric_area_obs1 = 0.0
		total_geometric_area_obs2 = 0.0		
		for line in open(self.param['foil_dimension_file']):
			cols = line.split() 
			if cols[0] != '#':
				foil_id                 = int(cols[0])
				foil_cone_center_radius = float(cols[7])
				foil_cone_angle_rad     = float(cols[8]) # rad
				foil_cone_angle_deg     = float(cols[9]) # degree
				foil_cone_top_radius    = float(cols[10])				
				foil_cone_bot_radius    = float(cols[11])				
				foil_cone_ann_wid_mm    = float(cols[12])				
				foil_geometric_area     = float(cols[18]) # mm2

				foil_parab_bot_radius   = float(cols[1]) # mm
				foil_parab_top_radius   = float(cols[5]) # mm 
				foil_parab_angle_rad    = float(cols[2])
				foil_parab_angle_deg    = scipy.degrees(foil_parab_angle_rad)

				foil_parab_bot_radius_cm = 0.1 * foil_parab_bot_radius
				foil_parab_top_radius_cm = 0.1 * foil_parab_top_radius
				self.foils.append(Foil(foil_id, 
					self.material, 
					foil_parab_top_radius_cm, foil_parab_bot_radius_cm, 
					foil_cone_angle_rad))
				#	foil_parab_angle_rad))
				self.foils[-1].roughness_A = self.param['roughness_A']
				self.foils[-1].n_spoke1    = self.param['spoke']['n_spoke1']
				self.foils[-1].n_spoke2    = self.param['spoke']['n_spoke2']
				self.foils[-1].spoke1_wid  = self.param['spoke']['spoke1_wid']
				self.foils[-1].spoke2_wid  = self.param['spoke']['spoke2_wid']
				#self.foils[-1].setUseFlag(self.use_foilid_list)
				self.foils[-1].setUseFlag(self.use_foilid_list)
				self.foils[-1].surface_model = self.param['surface_model']
				self.foils[-1].showProperty()

				total_geometric_area      += foil_geometric_area
				total_geometric_area_obs1 += foil_geometric_area * self.foils[-1].getSpokeObsculation()
				subtraction  = float(self.param['spoke']['n_spoke1']) * float(self.param['spoke']['spoke1_wid'])*10.0 * (self.foils[-1].top_radius_cm-self.foils[-1].bottom_radius_cm)*10.0
				subtraction += float(self.param['spoke']['n_spoke2']) * float(self.param['spoke']['spoke2_wid'])*10.0 * (self.foils[-1].top_radius_cm-self.foils[-1].bottom_radius_cm)*10.0
				value = foil_geometric_area - subtraction
				total_geometric_area_obs2 += value

		print "Total geometric area (cm2):",        total_geometric_area
		#print "after subtracting spoke obsculation: ", total_geometric_area_obs1
		print "after subtracting spoke obsculation: ", total_geometric_area_obs2
		#print self.n_spoke1, self.spoke1_wid
		#print self.n_spoke2, self.spoke2_wid
		#print self.foils[-1].top_radius_cm
		#print self.foils[0].bottom_radius_cm
		obst_1 = float(self.param['spoke']['n_spoke1']) * float(self.param['spoke']['spoke1_wid']) * (self.foils[-1].top_radius_cm-self.foils[0].bottom_radius_cm)
		obst_2 = float(self.param['spoke']['n_spoke2']) * float(self.param['spoke']['spoke2_wid']) * (self.foils[-1].top_radius_cm-self.foils[0].bottom_radius_cm)
		#print obst_1, obst_2
		#print total_geometric_area/100.0 - obst_1 - obst_2
		print "------------------------------------------------------------------------"

	def getFoil(self, foilid):
		return self.foils[foilid-1]

	def getEffectiveArea(self, energy_keV):
		Aeff_tot = 0.0
		for foil in self.foils:
			if foil.use_flag:
				#print "%d " % foil.foilid
				Aeff_tot += foil.getEffectiveArea(energy_keV)
		return scipy.real(Aeff_tot)
		#return Aeff_tot

	def generateEffectiveAreaFile(self):
		print("generateEffectiveAreaFile")

class Foil():
	def __init__(self, foilid, material, top_radius_cm, bottom_radius_cm, angle_rad):
		self.foilid   = foilid
		self.material = material

		self.top_radius_cm    = top_radius_cm
		self.bottom_radius_cm = bottom_radius_cm
		self.angle_rad = angle_rad
		self.angle_deg = scipy.degrees(self.angle_rad)
		self.geo_area_cm2 = self.getGeometricalArea()

		self.use_flag_mark = ' '			

		self.roughness_A = 0.0 
		self.n_spoke1    = 0
		self.n_spoke2    = 0
		self.spoke1_wid  = 0.0
		self.spoke2_wid  = 0.0		

	def isUseFoil(self, use_foilid_list):
			if self.foilid in use_foilid_list:
				return True
			else: 
				return False

	def setUseFlag(self, use_foilid_list):
			if self.isUseFoil(use_foilid_list):
				self.use_flag = True
				self.use_flag_mark = '*'
			else: 
				self.use_flag = False				
				self.use_flag_mark = ' '				

	def showProperty(self):
		print "Foil-ID %02d: %s %.7f %.7f %.7f %.7f %.8f   %s  (%s)" % (
			self.foilid,
			self.material.symbol,
			self.top_radius_cm,
			self.bottom_radius_cm,
			self.angle_rad,
			self.angle_deg,
			self.geo_area_cm2,
			self.use_flag_mark,
			self.surface_model)

	def getGeometricalArea(self):
		return pi*(self.top_radius_cm**2-self.bottom_radius_cm**2)

	"""
	def getSurfaceRoughness_Lalmbda_cm(self, lambda_cm):
		roughness_cm = self.roughness_A * A_cm
		numerator    = 4.0 * pi * roughness_cm * scipy.sin(self.angle_rad)
		denominator  = lambda_cm
		ratio = scipy.exp(-(numerator/denominator)**2)
		return ratio

	def getSurfaceRoughness(self, energy_keV):
		lambda_cm = convertEnergykeV2LambdaCM(energy_keV)
		return self.getSurfaceRoughness_Lalmbda_cm(lambda_cm)
	"""

	def getDebyeWallerFactor(self, energy_keV):
		return self.material.getDebyeWallerFactor(energy_keV, self.roughness_A, self.angle_deg)

	def getNevotCroceFactor(self, energy_keV):
		return self.material.getNevotCroceFactor(energy_keV, self.roughness_A, self.angle_deg)

	def getReflectivity(self, energy_keV):
		R  = self.material.getReflectivity(energy_keV, self.angle_rad)
		#ratio = self.getDebyeWallerFactor(energy_keV)
		if self.surface_model == "NC":
			ratio = self.getNevotCroceFactor(energy_keV)
		elif self.surface_model == "DW":
			ratio = self.getDebyeWallerFactor(energy_keV)
		return ratio * R

	def getSpokeObsculation(self):
		"""
		obsculation spoke width : cm (spoke1_wid, spoke2_wid)
		"""
		n_spoke1 = self.n_spoke1
		n_spoke2 = self.n_spoke2 
		spoke1_wid = self.spoke1_wid
		spoke2_wid = self.spoke2_wid 
		rmax = self.top_radius_cm
		rmin = self.bottom_radius_cm
		rwid = rmax - rmin
		area = self.geo_area_cm2
		obsculation = (area - n_spoke1*spoke1_wid*rwid - n_spoke2*spoke2_wid*rwid)/area
		return obsculation

	def getEffectiveArea(self, energy_keV):
		roughness_A = self.roughness_A	
		n_spoke1 = self.n_spoke1
		n_spoke2 = self.n_spoke2 
		spoke1_wid = self.spoke1_wid
		spoke2_wid = self.spoke2_wid 			
		obsculation = self.getSpokeObsculation()
		Reff        = self.getReflectivity(energy_keV)
		Ageo        = self.geo_area_cm2
		return obsculation * Ageo * Reff


class Material:
	def __init__(self, symbol, atomicWeight, density): 
		self.symbol       = symbol
		self.atomicWeight = atomicWeight
		self.density      = density 

		self.array_energy_eV  = scipy.array([],dtype='float')
		self.array_energy_keV = scipy.array([],dtype='float')	
		self.array_f1         = scipy.array([],dtype='float')
		self.array_f2         = scipy.array([],dtype='float')

		self.scatteringFactorFile = '%s/physics/sf/%s.nff' % (os.getenv('DATADIR_XRC'), self.symbol)
		self.readScatteringFactorFile(self.scatteringFactorFile)

	def readScatteringFactorFile(self, tablefile):
		flag_header = True
		for line in open(tablefile):
			if flag_header:
				flag_header = False
			else:
				cols = line.split()		
				self.array_energy_eV  = scipy.append(self.array_energy_eV, float(cols[0]))
				self.array_energy_keV = scipy.append(self.array_energy_keV, 0.001*float(cols[0]))
				self.array_f1         = scipy.append(self.array_f1, float(cols[1]))
				self.array_f2         = scipy.append(self.array_f2, float(cols[2]))
		self.interpolateScatteringFactor()

	def interpolateScatteringFactor(self):
		self.interpolate_keV_vs_f1 = interpolate.interp1d(
			self.array_energy_keV, self.array_f1,
			kind='linear', bounds_error=False, fill_value=0)
		self.interpolate_keV_vs_f2 = interpolate.interp1d(
			self.array_energy_keV, self.array_f2,
			kind='linear', bounds_error=False, fill_value=0)		
	
	def showProperty(self):
		print "---- Material -----"
		print "Material symbol : %s"   % self.symbol
		print "Atomic Weithg   : %.3f" % self.atomicWeight
		print "Density (g/cm^3): %.3f" % self.density
		print "-------------------"

	def getF1(self, energy_keV):
		return self.interpolate_keV_vs_f1(energy_keV)

	def getF2(self, energy_keV):
		return self.interpolate_keV_vs_f2(energy_keV)		

	def getOpticalConstant_Lalmbda_cm(self, f, lambda_cm):
		opticalConstant = ( Re_cm * N_A / 2.0 / pi ) * ( self.density / self.atomicWeight ) * lambda_cm * lambda_cm * f
		return opticalConstant 

	def getDelta_Lalmbda_cm(self, f1, lambda_cm):
		delta = self.getOpticalConstant_Lalmbda_cm(f1, lambda_cm)
		return delta

	def getBeta_Lalmbda_cm(self, f2, lambda_cm):
		beta = self.getOpticalConstant_Lalmbda_cm(f2, lambda_cm)
		return beta

	def getDelta(self, energy_keV):
		f1        = self.getF1(energy_keV)
		lambda_cm = convertEnergykeV2LambdaCM(energy_keV)
		return self.getDelta_Lalmbda_cm(f1, lambda_cm)

	def getBeta(self, energy_keV):
		f2        = self.getF2(energy_keV)
		lambda_cm = convertEnergykeV2LambdaCM(energy_keV)
		return self.getBeta_Lalmbda_cm(f2, lambda_cm)

	def getCriticalAngle_Delta(self, delta):
		"""
		output radian
		"""
		criticalAngle = scipy.sqrt(2.0*delta)
		return criticalAngle

	def getCriticalAngle(self, energy_keV):
		return self.getCriticalAngle_Delta(self.getDelta(energy_keV))

	def getReflectivity_Constant(self, beta, delta, theta, criticalAngle):
		A = beta/delta
		#if criticalAngle == 0:
		#	X = 0
		#else:
		#	X = theta/criticalAngle
		X = theta/criticalAngle
		H = X**2 + scipy.sqrt( (X**2-1)**2 + A**2 )
		numerator   = H - X * scipy.sqrt(2*(H-1.0))
		denominator = H + X * scipy.sqrt(2*(H-1.0))
		R = numerator / denominator
		return R

	def getReflectivity(self, energy_keV, theta):
		beta  = self.getBeta(energy_keV)
		delta = self.getDelta(energy_keV)
		criticalAngle = self.getCriticalAngle(energy_keV)
		return self.getReflectivity_Constant(beta, delta, theta, criticalAngle)

	"""
	def getDebyeWallerFactor(self, energy_keV, roughness_A, incident_angle_deg):
		lambda_A = convertEnergykeV2LambdaA(energy_keV)
		incident_angle_rad = scipy.radians(incident_angle_deg)
		q0_A   = 4.0 * pi * ReflectionIndexVacuum / lambda_A * scipy.sin(incident_angle_rad)
		factor = scipy.exp( - (q0_A * roughness_A)**2.0 )
		return factor
	"""

	def getDebyeWallerFactor(self, energy_keV, roughness_A, incident_angle_deg):
		roughness_cm = roughness_A * A_cm
		incident_angle_rad = scipy.radians(incident_angle_deg)
		numerator    = 4.0 * pi * roughness_cm * scipy.sin(incident_angle_rad)
		lambda_cm    = convertEnergykeV2LambdaCM(energy_keV)
		denominator  = lambda_cm
		factor = scipy.exp(-(numerator/denominator)**2)
		return factor	

	def getNevotCroceFactor(self, energy_keV, roughness_A, incident_angle_deg):
		delta = self.getDelta(energy_keV)		
		beta  = self.getBeta(energy_keV)		
		ReflectionIndex = (1-delta) + 1j * beta

		lambda_A           = convertEnergykeV2LambdaA(energy_keV)
		incident_angle_rad = scipy.radians(incident_angle_deg)
		q0_A   = 4.0 * pi * ReflectionIndexVacuum / lambda_A * scipy.sin(incident_angle_rad)
		q1_A   = 4.0 * pi / lambda_A * scipy.sqrt(ReflectionIndex**2.0-(scipy.cos(incident_angle_rad))**2.0)
		factor = scipy.exp( - q0_A * q1_A * roughness_A**2 )
		return scipy.absolute(factor)

	"""		
	def getNevotCroceFactor(self, energy_keV, roughness_A, incident_angle_deg):
		delta = self.getDelta(energy_keV)		
		beta  = self.getBeta(energy_keV)		
		ReflectionIndex = (1-delta) + 1j * beta

		lambda_A           = convertEnergykeV2LambdaA(energy_keV)
		incident_angle_rad = scipy.radians(incident_angle_deg)
		q0_A   = 4.0 * pi * ReflectionIndexVacuum / lambda_A * scipy.sin(incident_angle_rad)
		q1_A   = 4.0 * pi * ReflectionIndex / lambda_A * scipy.sin(incident_angle_rad)
		factor = scipy.exp( - q0_A * q1_A * roughness_A**2 )
		return scipy.absolute(factor)
		#return scipy.absolute(factor)
	"""

	"""
	def getNevotCroceFactor(self, energy_keV, roughness_A, incident_angle_deg):
		delta = self.getDelta(energy_keV)
		beta  = self.getBeta(energy_keV)
		reflectionIndex2 = (1-delta)**2.0 + beta**2.0 # square of reflection index 

		lambda_A           = convertEnergykeV2LambdaA(energy_keV)
		incident_angle_rad = scipy.radians(incident_angle_deg)
		q0_A   = 4.0 * pi * ReflectionIndexVacuum / lambda_A * scipy.sin(incident_angle_rad)
		q1_A   = 4.0 * pi * scipy.sqrt(reflectionIndex2-(scipy.cos(incident_angle_rad))**2.0)		
		factor = scipy.exp( - q0_A * q1_A * roughness_A**2.0 )
		out = scipy.sqrt(scipy.real(factor)**2+scipy.imag(factor)**2)
		return out
	"""


