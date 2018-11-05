#!/usr/bin/env python

"""
CLI for xspecresp.py 
"""

import click
import niresp.xspecresp

@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
	if ctx.invoked_subcommand is None:
		print(ctx.get_help())

@cli.command(help="Generate xspec-format rmf file.")
@click.argument("setup_yamlfile_path",type=click.Path(exists=True))
@click.argument("outrmffile",default='out/ni.rmf')
def generate_rmffile(setup_yamlfile_path,outrmffile):
	resp = niresp.xspecresp.XspecResponse(setup_yamlfile_path)	
	resp.generate_rmffile(outrmffile)

@cli.command(help="Generate xspec-format arf file.")
@click.argument("setup_yamlfile_path",type=click.Path(exists=True))
@click.argument("outarffile",default='out/ni.arf')
def generate_arffile(setup_yamlfile_path,outarffile):
	resp = niresp.xspecresp.XspecResponse(setup_yamlfile_path)	
	resp.generate_arffile(outarffile)

def main():
	cli()

if __name__ == "__main__":
	main()