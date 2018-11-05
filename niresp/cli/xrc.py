#!/usr/bin/env python

"""
CLI for xrc.py 
"""

import click
import numpy as np
import niresp.xrc

@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
	if ctx.invoked_subcommand is None:
		print(ctx.get_help())

@cli.command(help="XRC.")
@click.argument("setup_yamlfile_path",type=click.Path(exists=True),default='data/xrc/xrc.yaml')
def plotEffectiveArea(setup_yamlfile_path):
	print("plot effective area")
	#xrc = niresp.xrc.XrayConcentrator(setup_yamlfile_path)	
	#energy_keV_array = np.arange(0.1,15,0.1)
	#effective_area = [xrc.getEffectiveArea(energy_keV) for energy_keV in energy_keV_array]
	#print(effective_area)

@cli.command(help="get effective area of NICER with an input yaml file.")
@click.argument("energy_kev")
@click.argument("setup_yamlfile_path",type=click.Path(exists=True),default='data/xrc/xrc.yaml')
def getEffectiveArea(energy_kev,setup_yamlfile_path):
	xrc = niresp.xrc.XrayConcentrator(setup_yamlfile_path)	
	print("Total effective area {:.4f} cm2 at {} keV".format(xrc.getEffectiveArea(float(energy_kev)),energy_kev))

@cli.command(help="XRC.")
@click.argument("setup_yamlfile_path",type=click.Path(exists=True),default='data/xrc/xrc.yaml')
def plotReflectivity(setup_yamlfile_path):
	print("plot reflectivity")

def main():
	cli()

if __name__ == "__main__":
	main()