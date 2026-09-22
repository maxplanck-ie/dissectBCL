import os

import rich_click as click
from rich import print
from rich.table import Table

from dissectBCL.misc import getConf, getVersion
from wd40.fex import fex as fex_upload
from wd40.release import rel as release
from wd40.reset import reset as reset_outLane

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
            },
        }
    ]
}

click.rich_click.COMMAND_GROUPS = {
    "wd40": [
        {
            "name": "Main commands",
            "commands": ["rel", "reset", "fex", "help"],
        }
    ]
}

COMMAND_HELP = {
    "rel": (
        "wd40 rel [flowcell]",
        "Release a finished flowcell to periphery storage: chmod/chgrp the "
        "flowcell, project, FASTQC, and Analysis folders, and push filepaths "
        "to Parkour2. Run after BigRedButton has set analysis.done.",
    ),
    "reset": (
        "wd40 reset [outLane]",
        "Strip an outLane dir under /rapidus back to just its "
        "SampleSheet/RunManifest, deleting demux output and done-flags. Use "
        "it to hand-edit the samplesheet (index mask, mismatches, I5/dual vs "
        "single index) and redemux, without re-copying from the flowcell's "
        "read-only source directory.",
    ),
    "fex": (
        "wd40 fex [project]",
        "Upload a dissectBCL project to FEX as an RO-Crate archive. Fetches "
        "comprehensive ISA-profile metadata from parkour-test (latest fixes), "
        "enriches it with FASTQ file entities and md5 checksums, and streams "
        "the zip to fexsend without writing to disk. Project name must match "
        "Project_XXXX_User_PI format.",
    ),
    "help": (
        "wd40 help",
        "Show this list of subcommands and when to reach for each.",
    ),
}


@click.group()
@click.option(
    "--configpath",
    show_default=True,
    required=False,
    default=os.path.expanduser("~/configs/dissectBCL_prod.ini"),
    help="config file location",
    type=click.Path(exists=True),
)
@click.option(
    "--debug/--no-debug",
    "-d/-n",
    default=False,
    show_default=True,
    help="Show the debug log messages",
)
@click.version_option(getVersion("dissectBCL"), prog_name="wd40")
@click.pass_context
def cli(ctx, configpath, debug):
    ctx.ensure_object(dict)
    ctx.obj["DEBUG"] = debug
    ctx.obj["configpath"] = configpath
    # populate ctx from config.
    # For release:
    cnf = getConf(configpath, quickload=True)
    ctx.obj["prefixDir"] = cnf["Dirs"]["piDir"]
    ctx.obj["piList"] = cnf["Internals"]["PIs"]
    ctx.obj["postfixDir"] = cnf["Internals"]["seqDir"]
    #    ctx.obj['solDir'] = cnf['Dirs']['baseDir']
    ctx.obj["parkourURL"] = cnf["parkour"]["URL"]
    ctx.obj["parkourAuth"] = (cnf["parkour"]["user"], cnf["parkour"]["password"])
    ctx.obj["parkourCert"] = cnf["parkour"]["cert"]
    ctx.obj["fexBool"] = cnf["Internals"].getboolean("fex")
    ctx.obj["fromAddress"] = cnf["communication"]["fromAddress"]


@cli.command()
@click.argument("flowcell", default="./", type=click.Path(exists=True))
@click.pass_context
def rel(ctx, flowcell):
    """Releases a flowcell."""
    release(
        flowcell,
        ctx.obj["piList"],
        ctx.obj["prefixDir"],
        ctx.obj["postfixDir"],
        ctx.obj["parkourURL"],
        ctx.obj["parkourAuth"],
        ctx.obj["parkourCert"],
        ctx.obj["fexBool"],
        ctx.obj["fromAddress"],
    )


@cli.command()
@click.argument("outlane", default="./", type=click.Path(exists=True))
def reset(outlane):
    """Strips an outLane dir back to its SampleSheet/RunManifest, for hand-editing."""
    reset_outLane(outlane)


@cli.command()
@click.argument("project", type=click.Path(exists=True))
@click.option(
    "--parkour-url",
    default=None,
    help="Override parkour URL (default: parkour-test for latest fixes)",
)
@click.pass_context
def fex(ctx, project, parkour_url):
    """Upload a project to FEX as an RO-Crate archive."""
    config = getConf(ctx.obj["configpath"], quickload=True)
    fex_upload(project, config, ctx.obj["fromAddress"], parkour_url)


@cli.command(name="help")
def help_cmd():
    """Lists all subcommands and when to use them."""
    table = Table(title="wd40 subcommands")
    table.add_column("Usage", style="cyan")
    table.add_column("When to use it")
    for usage, when in COMMAND_HELP.values():
        table.add_row(usage, when)
    print(table)
