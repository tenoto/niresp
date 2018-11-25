#!/usr/bin/env python

"""
CLI for xrc.py 
"""

import click
import niresp.cubesat

@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
	if ctx.invoked_subcommand is None:
		print(ctx.get_help())

@cli.command(help="Plot and write effective area curve for a Cubesat.")
@click.argument("setup_yamlfile_path",type=click.Path(exists=True),
	default='data/cubesat/cubesat.yaml')
@click.option('--foil_dimension_file', '-f', 
	default='data/cubesat/foilparameter/v181108/foilpara/foilpara_fl30cm.txt')
@click.option('--maximum_incident_angle_deg', '-a', default=3.0)
@click.option('--maximum_diameter_cm', '-d', default=9.0)
@click.option('--scaling_factor', '-d', default=1.0)
@click.option('--outputfile', '-o', default='cubesat_effectivearea.txt')
def generateEffectiveAreaFile(setup_yamlfile_path,
	foil_dimension_file,
	maximum_incident_angle_deg,maximum_diameter_cm,
	scaling_factor,
	outputfile):
	cubesat_xrc = niresp.cubesat.CubeSatXRC(setup_yamlfile_path)
	cubesat_xrc.set_foil_dimension_file(foil_dimension_file)
	cubesat_xrc.set_maximum_incident_angle_deg(maximum_incident_angle_deg)
	cubesat_xrc.set_maximum_diameter_cm(maximum_diameter_cm)	
	cubesat_xrc.set_scaling_factor(scaling_factor)		
	cubesat_xrc.setMaterial()		
	cubesat_xrc.showProperty()
	cubesat_xrc.setFoilDimension()
	cubesat_xrc.generateEffectiveAreaFile(outfile=outputfile)

@cli.command(help="Calculate an effective area at the input energy.")
@click.argument("setup_yamlfile_path",type=click.Path(exists=True),
	default='data/cubesat/cubesat.yaml')
@click.argument("energy_kev",type=float, nargs=1)
def getEffectiveArea(setup_yamlfile_path,energy_kev):
	print("get effective area")
	cubesat_xrc = niresp.cubesat.CubeSatXRC(setup_yamlfile_path)
	print(cubesat_xrc.getEffectiveArea(energy_kev))

@cli.command(help="Generate rmf file for a Cubesat.")
@click.argument("setup_yamlfile_path",type=click.Path(exists=True),
	default='data/cubesat/cubesat.yaml')
@click.argument("outrmffile",type=str, nargs=1,default='gwcubesat.rmf')
def generateRmffile(setup_yamlfile_path,outrmffile):
	print("gnerate arf file for a CubeSat.")
	cubesat = niresp.cubesat.CubeSatXspecResponse(setup_yamlfile_path)
	cubesat.generate_rmffile(outrmffile)

@cli.command(help="Generate arf file for a Cubesat.")
@click.argument("setup_yamlfile_path",type=click.Path(exists=True),
	default='data/cubesat/cubesat.yaml')
@click.argument("inareafile",type=click.Path(exists=True))
@click.argument("outarffile",type=str, nargs=1,default='gwcubesat.arf')
def generateArffile(setup_yamlfile_path,inareafile,outarffile):
	print("gnerate arf file for a CubeSat.")
	cubesat = niresp.cubesat.CubeSatXspecResponse(setup_yamlfile_path)
	cubesat.generate_arffile(inareafile,outarffile)

def main():
	cli()

if __name__ == "__main__":
	main()