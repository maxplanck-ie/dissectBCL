import os

import rich_click as click
from rich import print

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
            "options": ["--configpath", "--help", "--version", "--debug"],
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


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "--configpath",
    show_default=True,
    required=False,
    default=os.path.expanduser("~/configs/dissectBCL_prod.ini"),
    help="Path to the dissectBCL config .ini file",
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
    """wd40: post-demux helpers (release, reset, FEX upload) for dissectBCL runs."""
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
    """Release a finished flowcell to periphery storage.

    chmod/chgrp the flowcell, project, FASTQC, and Analysis folders, and push
    filepaths to Parkour2. Run after BigRedButton has set analysis.done.

    FLOWCELL: path to the flowcell directory (default: current directory).
    Must contain analysis.done, set by BigRedButton.
    """
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
    """Strip an outLane dir back to its SampleSheet/RunManifest, for hand-editing.

    Strips an outLane dir under /rapidus back to just its SampleSheet/RunManifest,
    deleting demux output and done-flags. Use it to hand-edit the samplesheet
    (index mask, mismatches, I5/dual vs single index) and redemux, without
    re-copying from the flowcell's read-only source directory.

    OUTLANE: path to the outLane directory (default: current directory).
    """
    reset_outLane(outlane)


@cli.command()
@click.argument("project", type=click.Path(exists=True))
@click.option(
    "--parkour-url",
    default=None,
    help="Override parkour URL (default: parkour.URL from the config file)",
)
@click.pass_context
def fex(ctx, project, parkour_url):
    """Upload a project to FEX as an RO-Crate archive.

    Fetches comprehensive ISA-profile metadata from Parkour, enriches it with
    FASTQ file entities and md5 checksums, and streams the zip to fexsend.
    Projects >= 4 GiB are zipped to a temp file first, removed afterwards.

    PROJECT: path to the project directory, named Project_XXXX_User_PI.
    """
    config = getConf(ctx.obj["configpath"], quickload=True)
    fex_upload(project, config, ctx.obj["fromAddress"], parkour_url)


@cli.command(name="help")
@click.pass_context
def help_cmd(ctx):
    """Show the full help page of every subcommand, as if running each with -h."""
    for name, cmd in cli.commands.items():
        if name == "help":
            continue
        sub = cmd.context_class(cmd, info_name=name, parent=ctx.parent)
        click.echo(sub.get_help())
