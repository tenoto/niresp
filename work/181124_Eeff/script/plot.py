#!/usr/bin/env python

import pandas 
import astropy.io.fits as pyfits

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages 
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib import rc
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = '14'
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['xtick.top'] = 'True'
mpl.rcParams['ytick.right'] = 'True'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rcParams['axes.grid'] = 'True'
mpl.rcParams['axes.xmargin'] = '.05' #'.05'
mpl.rcParams['axes.ymargin'] = '.05'
mpl.rcParams['savefig.facecolor'] = 'None'
mpl.rcParams['savefig.edgecolor'] = 'None'
mpl.rcParams['savefig.bbox'] = 'tight'

class gwCubeSatArffile():
	def __init__(self,arffile):
		self.arffile = arffile

		self.focal_length_cm = (float(self.arffile.split('_fl')[-1].split('cm')[0].replace('p','.')))

		self.hdu = pyfits.open(arffile)
		self.energy_keV  = self.hdu['SPECRESP'].data['ENERGY']
		self.area_cm2 = self.hdu['SPECRESP'].data['SPECRESP']
		self.num_of_foil = self.hdu['SPECRESP'].header['NUMFOILS']
		self.scaleing_factor = self.hdu['SPECRESP'].header['SCALINGF']
		self.label = 'FL=%.1f cm (%d foils)' % (self.focal_length_cm,self.num_of_foil)

class nicerArffile():
	def __init__(self,arffile):
		self.arffile = arffile

		self.hdu = pyfits.open(arffile)
		self.energy_keV  = self.hdu['SPECRESP'].data['ENERGY']
		self.area_cm2 = self.hdu['SPECRESP'].data['SPECRESP']

# ===================================
# Main 
# ===================================
pp = PdfPages('work/181124_Eeff/effective_area/gwcubesat_xrc_enoto_v181124a.pdf')


arf_fl24p9_i3p0_s1p00 = gwCubeSatArffile('work/181124_Eeff/effective_area/foilpara_fl24p92cm_i3p0_s1p00_v181124a.arf')
arf_fl50p0_i3p0_s1p00 = gwCubeSatArffile('work/181124_Eeff/effective_area/foilpara_fl50cm_i3p0_s1p00_v181124a.arf')
arf_fl24p9_i3p0_s0p85 = gwCubeSatArffile('work/181124_Eeff/effective_area/foilpara_fl24p92cm_i3p0_s0p85_v181124a.arf')
arf_fl50p0_i3p0_s0p85 = gwCubeSatArffile('work/181124_Eeff/effective_area/foilpara_fl50cm_i3p0_s0p85_v181124a.arf')

nicer_00 = nicerArffile('/Users/enoto/Dropbox/enoto/research/nicer/response/180305_ver1.02/archive/nicer_resp_ver1.02_detid/nicer_resp_ver1.02_detid/ni_detid00_onaxis_v1.02.arf')

list_s1p00 = [arf_fl24p9_i3p0_s1p00,arf_fl50p0_i3p0_s1p00]
list_s0p85 = [arf_fl24p9_i3p0_s0p85,arf_fl50p0_i3p0_s0p85]
list_color = ['red','blue']

fig = plt.figure()
for i in range(len(list_s1p00)):
	plt.plot(list_s1p00[i].energy_keV,list_s1p00[i].area_cm2,label=list_s1p00[i].label,	
		color=list_color[i])
	plt.plot(list_s0p85[i].energy_keV,list_s0p85[i].area_cm2,color=list_color[i],lw=1,ls='--')
plt.plot(nicer_00.energy_keV,nicer_00.area_cm2,ls='-.',lw=1,color='#AEB6BF',label='NICER single module')	
plt.xlabel('Energy (keV)')
plt.ylabel('Effective Area (cm$^2$)')
plt.xscale('log')
plt.legend(title='Incident angle < 3.0 deg',fontsize=10,title_fontsize='small',loc='upper left')
plt.ylim(0,35)
plt.xlim(0.15,10)
plt.title('XRC effective area with filter attenuations (GWCubeSat, v181124a)',fontsize='small')
ax = plt.gca()
ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
pp.savefig(fig)

df_rate = pandas.read_table('work/181124_Eeff/effective_area/fl_vs_rate_crab_scox1.txt',
	header=0,index_col=None,skiprows=[0,1,2,3],delim_whitespace=True)
print(df_rate)
fig = plt.figure()
plt.plot(df_rate['fl'],df_rate['crab_rate'],'o',label='Crab')
plt.plot(df_rate['fl'],df_rate['scox1_rate'],'s',label='Sco X-1')
plt.xlabel('Focal length (cm)')
plt.ylabel('Count rate (cps)')
plt.legend(loc='upper left',fontsize=10,title_fontsize='small')
plt.ylim(0,300)
plt.xlim(0,70)
plt.title('Count rate in the 0.2-3.0 keV band (GWCubeSat, v181124a)',fontsize='small')
#ax = plt.gca()
#ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
pp.savefig(fig)




"""
# ===================================
# Page 1 
# ===================================
df_fl108p5cm_1p4= pandas.read_table('work/181124_Eeff/effective_area/foilpara_fl108p5cm_i1p4_s1p00_v181124a.txt',
	header=None,names=('energy_keV','area_cm2'),index_col=None,sep=' ',skiprows=1)
df_fl60cm_1p4 = pandas.read_table('work/181124_Eeff/effective_area/foilpara_fl60cm_i1p4_s1p00_v181124a.txt',
	header=None,names=('energy_keV','area_cm2'),index_col=None,sep=' ',skiprows=1)
df_fl50cm_1p4 = pandas.read_table('work/181124_Eeff/effective_area/foilpara_fl50cm_i1p4_s1p00_v181124a.txt',
	header=None,names=('energy_keV','area_cm2'),index_col=None,sep=' ',skiprows=1)
df_fl40cm_1p4 = pandas.read_table('work/181124_Eeff/effective_area/foilpara_fl40cm_i1p4_s1p00_v181124a.txt',
	header=None,names=('energy_keV','area_cm2'),index_col=None,sep=' ',skiprows=1)
df_fl40cm_1p4 = pandas.read_table('work/181124_Eeff/effective_area/foilpara_fl30cm_i1p4_s1p00_v181124a.txt',
	header=None,names=('energy_keV','area_cm2'),index_col=None,sep=' ',skiprows=1)

fig = plt.figure()
plt.plot(df_fl60cm_1p4['energy_keV'],df_fl60cm_1p4['area_cm2'],label='FL60 cm (7 foils)')
plt.plot(df_fl50cm_1p4['energy_keV'],df_fl50cm_1p4['area_cm2'],label='FL50 cm (5 foils)')
plt.plot(df_fl40cm_1p4['energy_keV'],df_fl40cm_1p4['area_cm2'],label='FL40 cm (2 foils)')
plt.plot(df_fl108p5cm_1p4['energy_keV'],df_fl108p5cm_1p4['area_cm2'],label='FL108.5 cm (24 foils)')
plt.xlabel('Energy (keV)')
plt.ylabel('Effective Area (cm$^2$)')
plt.xscale('log')
plt.legend(title='Incident angle < 1.4 deg',fontsize=10,title_fontsize='small')
plt.ylim(0,120)
plt.xlim(0.1,20)
plt.title('XRC area curve without any filter attenuation (GWCubeSat, v181124a)',fontsize='small')
ax = plt.gca()
ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
pp.savefig(fig)

# ===================================
# Page 2
# ===================================
df_fl108p5cm_3p0 = pandas.read_table('work/181124_Eeff/effective_area/foilpara_fl108p5cm_3p0_v181124a.txt',
	header=None,names=('energy_keV','area_cm2'),index_col=None,sep=' ',skiprows=1)
df_fl60cm_3p0 = pandas.read_table('work/181124_Eeff/effective_area/foilpara_fl60cm_3p0_v181124a.txt',
	header=None,names=('energy_keV','area_cm2'),index_col=None,sep=' ',skiprows=1)
df_fl50cm_3p0 = pandas.read_table('work/181124_Eeff/effective_area/foilpara_fl50cm_3p0_v181124a.txt',
	header=None,names=('energy_keV','area_cm2'),index_col=None,sep=' ',skiprows=1)
df_fl40cm_3p0 = pandas.read_table('work/181124_Eeff/effective_area/foilpara_fl40cm_3p0_v181124a.txt',
	header=None,names=('energy_keV','area_cm2'),index_col=None,sep=' ',skiprows=1)
df_fl30cm_3p0 = pandas.read_table('work/181124_Eeff/effective_area/foilpara_fl30cm_3p0_v181124a.txt',
	header=None,names=('energy_keV','area_cm2'),index_col=None,sep=' ',skiprows=1)

fig = plt.figure()
plt.plot(df_fl60cm_3p0['energy_keV'],df_fl60cm_3p0['area_cm2'],label='FL60 cm (16 foils)')
plt.plot(df_fl50cm_3p0['energy_keV'],df_fl50cm_3p0['area_cm2'],label='FL50 cm (12 foils)')
plt.plot(df_fl40cm_3p0['energy_keV'],df_fl40cm_3p0['area_cm2'],label='FL40 cm (8 foils)')
plt.plot(df_fl30cm_3p0['energy_keV'],df_fl30cm_3p0['area_cm2'],label='FL30 cm (4 foils)')
#plt.plot(df_fl108p5cm_3p0['energy_keV'],df_fl108p5cm_3p0['area_cm2'],label='FL108.5 cm (24 foils)')
plt.xlabel('Energy (keV)')
plt.ylabel('Effective Area (cm$^2$)')
plt.xscale('log')
plt.legend(title='Incident angle < 3 deg',fontsize=10,title_fontsize='small')
plt.ylim(0,120)
plt.xlim(0.1,20)
plt.title('XRC area curve without any filter attenuation (GWCubeSat, v181124a)',fontsize='small')
ax = plt.gca()
ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
pp.savefig(fig)

# ===================================
# Page 3
# ===================================
hdu_fl60cm_1p4 = pyfits.open('work/181124_Eeff/effective_area/foilpara_fl60cm_1p4_v181124a.arf')
hdu_fl50cm_1p4 = pyfits.open('work/181124_Eeff/effective_area/foilpara_fl50cm_1p4_v181124a.arf')
hdu_fl40cm_1p4 = pyfits.open('work/181124_Eeff/effective_area/foilpara_fl40cm_1p4_v181124a.arf')
hdu_fl108p5cm_1p4 = pyfits.open('work/181124_Eeff/effective_area/foilpara_fl108p5cm_1p4_v181124a.arf')

fig = plt.figure()
plt.plot(
	hdu_fl60cm_1p4['SPECRESP'].data['ENERGY'],
	hdu_fl60cm_1p4['SPECRESP'].data['SPECRESP'],
	label='FL60 cm (7 foils)')
plt.plot(
	hdu_fl50cm_1p4['SPECRESP'].data['ENERGY'],
	hdu_fl50cm_1p4['SPECRESP'].data['SPECRESP'],
	label='FL50 cm (5 foils)')
plt.plot(
	hdu_fl40cm_1p4['SPECRESP'].data['ENERGY'],
	hdu_fl40cm_1p4['SPECRESP'].data['SPECRESP'],
	label='FL40 cm (2 foils)')
plt.plot(
	hdu_fl108p5cm_1p4['SPECRESP'].data['ENERGY'],
	hdu_fl108p5cm_1p4['SPECRESP'].data['SPECRESP'],
	label='FL108.5 cm (24 foils)')
plt.xlabel('Energy (keV)')
plt.ylabel('Effective Area (cm$^2$)')
plt.xscale('log')
plt.legend(title='Incident angle < 1.4 deg',fontsize=10,title_fontsize='small')
plt.ylim(0,50)
plt.xlim(0.1,20)
plt.title('XRC area curve witt filter attenuations (GWCubeSat, v181124a)',fontsize='small')
ax = plt.gca()
ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
pp.savefig(fig)

# ===================================
# Page 4
# ===================================
hdu_fl60cm_3p0 = pyfits.open('work/181124_Eeff/effective_area/foilpara_fl60cm_3p0_v181124a.arf')
hdu_fl50cm_3p0 = pyfits.open('work/181124_Eeff/effective_area/foilpara_fl50cm_3p0_v181124a.arf')
hdu_fl40cm_3p0 = pyfits.open('work/181124_Eeff/effective_area/foilpara_fl40cm_3p0_v181124a.arf')
hdu_fl30cm_3p0 = pyfits.open('work/181124_Eeff/effective_area/foilpara_fl30cm_3p0_v181124a.arf')
hdu_fl20cm_3p0 = pyfits.open('work/181124_Eeff/effective_area/foilpara_fl20cm_3p0_v181124a.arf')
hdu_fl108p5cm_3p0 = pyfits.open('work/181124_Eeff/effective_area/foilpara_fl108p5cm_3p0_v181124a.arf')

fig = plt.figure()
plt.plot(
	hdu_fl60cm_3p0['SPECRESP'].data['ENERGY'],
	hdu_fl60cm_3p0['SPECRESP'].data['SPECRESP'],
	label='FL60 cm (16 foils)')
plt.plot(
	hdu_fl50cm_3p0['SPECRESP'].data['ENERGY'],
	hdu_fl50cm_3p0['SPECRESP'].data['SPECRESP'],
	label='FL50 cm (12 foils)')
plt.plot(
	hdu_fl40cm_3p0['SPECRESP'].data['ENERGY'],
	hdu_fl40cm_3p0['SPECRESP'].data['SPECRESP'],
	label='FL40 cm (8 foils)')
plt.plot(
	hdu_fl30cm_3p0['SPECRESP'].data['ENERGY'],
	hdu_fl30cm_3p0['SPECRESP'].data['SPECRESP'],
	label='FL30 cm (4 foils)')
#plt.plot(
#	hdu_fl108p5cm_3p0['SPECRESP'].data['ENERGY'],
#	hdu_fl108p5cm_3p0['SPECRESP'].data['SPECRESP'],
#	label='FL108.5 cm (24 foils)')
plt.xlabel('Energy (keV)')
plt.ylabel('Effective Area (cm$^2$)')
plt.xscale('log')
plt.legend(title='Incident angle < 3 deg',fontsize=10,title_fontsize='small')
plt.ylim(0,50)
plt.xlim(0.1,20)
plt.title('XRC area curve witt filter attenuations (GWCubeSat, v181124a)',fontsize='small')
ax = plt.gca()
ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
pp.savefig(fig)

# ===================================
# Page 5
# ===================================
row = 210
fl_3p0_array = []
fl_1p4_array = []
area_cm2_3p0_array = []
area_cm2_1p4_array = []
print(hdu_fl30cm_3p0['SPECRESP'].data['ENERGY'][row])
fl_3p0_array.append(20.0);area_cm2_3p0_array.append(hdu_fl20cm_3p0['SPECRESP'].data['SPECRESP'][row])
fl_3p0_array.append(30.0);area_cm2_3p0_array.append(hdu_fl30cm_3p0['SPECRESP'].data['SPECRESP'][row])
fl_3p0_array.append(40.0);area_cm2_3p0_array.append(hdu_fl40cm_3p0['SPECRESP'].data['SPECRESP'][row])
fl_3p0_array.append(50.0);area_cm2_3p0_array.append(hdu_fl50cm_3p0['SPECRESP'].data['SPECRESP'][row])
fl_3p0_array.append(60.0);area_cm2_3p0_array.append(hdu_fl60cm_3p0['SPECRESP'].data['SPECRESP'][row])

fl_1p4_array.append(20.0);area_cm2_1p4_array.append(0.0)
fl_1p4_array.append(30.0);area_cm2_1p4_array.append(0.0)
fl_1p4_array.append(40.0);area_cm2_1p4_array.append(hdu_fl40cm_1p4['SPECRESP'].data['SPECRESP'][row])
fl_1p4_array.append(50.0);area_cm2_1p4_array.append(hdu_fl50cm_1p4['SPECRESP'].data['SPECRESP'][row])
fl_1p4_array.append(60.0);area_cm2_1p4_array.append(hdu_fl60cm_1p4['SPECRESP'].data['SPECRESP'][row])
fl_1p4_array.append(108.5);area_cm2_1p4_array.append(hdu_fl108p5cm_1p4['SPECRESP'].data['SPECRESP'][row])

fig = plt.figure()
plt.plot(fl_3p0_array,area_cm2_3p0_array,'o',label='Incident angle < 3 deg')
plt.plot(fl_1p4_array,area_cm2_1p4_array,'s',label='Incident angle < 1.4 deg')
plt.xlabel('Focal length (cm)')
plt.ylabel('Effective area (cm$^{2}$)')
plt.legend(loc='upper left',fontsize=10,title_fontsize='small')
plt.ylim(0,50)
plt.xlim(0,110)
plt.title('Effective area at 1.15 keV after filter attenuations (GWCubeSat, v181124a)',fontsize='small')
#ax = plt.gca()
#ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
pp.savefig(fig)

# ===================================
# Page 6
# ===================================
df_rate = pandas.read_table('work/181124_Eeff/data/crab_rate.txt',
	header=0,index_col=None,skiprows=[1,2],delim_whitespace=True)
print(df_rate)
fig = plt.figure()
plt.plot(df_rate['FL'],df_rate['rate_3p0'],'o',label='Incident angle < 3 deg')
plt.plot(df_rate['FL'],df_rate['rate_1p4'],'s',label='Incident angle < 1.4 deg')
plt.xlabel('Focal length (cm)')
plt.ylabel('Count rate of the Crab (cps)')
plt.legend(loc='upper left',fontsize=10,title_fontsize='small')
plt.ylim(0,250)
plt.xlim(0,110)
plt.title('Crab count rate in the 0.5-3.0 keV band (GWCubeSat, v181124a)',fontsize='small')
#ax = plt.gca()
#ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
pp.savefig(fig)
"""

pp.close()