import rich_click as click
from rich import print
from importlib.metadata import version
import os
import sys
from dissectBCL.misc import getConf
from wd40.release import rel as release
from wd40.cat import catRun
from wd40.diagnose import diagnose

can_string = "[red]            ___ \n[/red]"
can_string += "[red]           |___|--------[/red]\n"
can_string += "[blue]           |   |[/blue]\n"
can_string += "[blue]           | [yellow]W[/yellow] |[/blue]\n"
can_string += "[blue]    WD40   | [yellow]D[/yellow] |  Kriechöl[/blue]\n"
can_string += "[blue]           | [yellow]4[/yellow] |[/blue]\n"
can_string += "[blue]           | [yellow]0[/yellow] |[/blue]\n"
can_string += "[blue]           |___|[/blue]\n"
print(can_string)

click.rich_click.OPTION_GROUPS = {
    "wd40": [
        {
            "name": "Options",
            "Options": ["--configpath", "--help", "--version", "--debug"],
            "table_styles": {
                "row_styles": ["cyan", "cyan", "cyan", "cyan"],
            }
        }
    ]
}

click.rich_click.COMMAND_GROUPS = {
    "wd40": [
        {
            "name": "Main commands",
            "commands": ["rel", "cat", "diag"],
        }
    ]
}


@click.group()
@click.option(
    "--configpath",
    show_default=True,
    required=False,
    default=os.path.expanduser('~/configs/dissectBCL_prod.ini'),
    help="config file location",
    type=click.Path(exists=True)
)
@click.option(
    "--debug/--no-debug",
    "-d/-n",
    default=False,
    show_default=True,
    help="Show the debug log messages",
)
@click.version_option(
    version("dissectBCL"),
    prog_name='wd40'
)
@click.pass_context
def cli(ctx, configpath, debug):
    ctx.ensure_object(dict)
    ctx.obj['DEBUG'] = debug
    ctx.obj['configpath'] = configpath
    # populate ctx from config.
    # For release:
    cnf = getConf(configpath)
    ctx.obj['prefixDir'] = cnf['Dirs']['piDir']
    ctx.obj['piList'] = cnf['Internals']['PIs']
    ctx.obj['postfixDir'] = cnf['Internals']['seqDir']
    ctx.obj['fastqDir'] = cnf['Dirs']['outputDir']
    ctx.obj['solDir'] = cnf['Dirs']['baseDir']


@cli.command()
@click.argument(
    "flowcell",
    default='./',
    type=click.Path(exists=True)
)
@click.pass_context
def rel(ctx, flowcell):
    """Releases a flowcell."""
    release(
        flowcell,
        ctx.obj['piList'],
        ctx.obj['prefixDir'],
        ctx.obj['postfixDir']
    )


@cli.command()
@click.option(
    "--flowcells",
    "-f",
    required=True,
    help='Specify 2 flowcell folders, e.g. -f flowcell1 -f flowcell2',
    multiple=True
)
@click.option(
    "--project",
    "-p",
    required=True,
    help='project folder. only 1 allowed.'
)
@click.option(
    "--output",
    "-o",
    required=True,
    help='folder to write output into.'
)
@click.pass_context
def cat(ctx, flowcells, project, output):
    """combine fastq files of a project sequenced on multiple flow cells."""
    if len(flowcells) != 2:
        sys.exit("Please specify exactly two flowcells..")
    project = os.path.basename(project)
    p1 = os.path.join(
        ctx.obj['fastqDir'],
        os.path.basename(flowcells[0]),
        project
    )
    p2 = os.path.join(
        ctx.obj['fastqDir'],
        os.path.basename(flowcells[1]),
        project
    )
    if not os.path.exists(p1):
        sys.exit('{} not found.'.format(p1))
    if not os.path.exists(p2):
        sys.exit('{} not found.'.format(p2))
    catRun(project, p1, p2, os.path.abspath(output))


@cli.command()
@click.argument(
    "flowcell",
    default='./',
    type=click.Path(exists=True)
)
@click.pass_context
def diag(ctx, flowcell):
    """Diagnose a flowcell."""
    diagnose(flowcell, ctx.obj['solDir'])
