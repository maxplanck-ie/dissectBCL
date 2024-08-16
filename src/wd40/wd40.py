import rich_click as click
from rich import print
from importlib.metadata import version
import os
from dissectBCL.misc import getConf
from wd40.release import rel as release

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
            "commands": ["rel"],
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
    cnf = getConf(configpath, quickload=True)
    ctx.obj['prefixDir'] = cnf['Dirs']['piDir']
    ctx.obj['piList'] = cnf['Internals']['PIs']
    ctx.obj['postfixDir'] = cnf['Internals']['seqDir']
    ctx.obj['fastqDir'] = cnf['Dirs']['outputDir']
    ctx.obj['solDir'] = cnf['Dirs']['baseDir']
    ctx.obj['parkourURL'] = cnf['parkour']['URL']
    ctx.obj['parkourAuth'] = (
                            cnf['parkour']['user'],
                            cnf['parkour']['password']
                            )
    ctx.obj['parkourCert'] = cnf['parkour']['cert']
    ctx.obj['fexBool'] = cnf['Internals'].getboolean('fex')
    ctx.obj['fromAddress'] = cnf['communication']['fromAddress']


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
        ctx.obj['postfixDir'],
        ctx.obj['parkourURL'],
        ctx.obj['parkourAuth'],
        ctx.obj['parkourCert'],
        ctx.obj['fexBool'],
        ctx.obj['fromAddress']
    )
