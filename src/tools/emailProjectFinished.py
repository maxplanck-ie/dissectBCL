#!/package/anaconda3/bin/python3
import argparse
import glob
import os
import re
import smtplib
import subprocess as sp
import sys
from email.mime.text import MIMEText
from urllib.parse import quote, unquote, urlsplit, urlunsplit

import requests

from dissectBCL.misc import getConf, projectPI


def getContactDetails(projectID, config):
    """
    Retrieve user data from a given sequencing request
    """
    res = requests.get(
        config["parkour"]["URL"]
        + "/api/requests/"
        + projectID
        + "/get_contact_details",
        auth=(config["parkour"]["user"], config["parkour"]["password"]),
        verify=config["parkour"]["cert"],
    )
    if res.status_code != 200:
        raise RuntimeError(f"API error: {res.json()}")
    return res.json()


def getProjectIDs(projects, config, forcePI=None):
    IDs = []
    for p in projects:
        # Sanity check
        assert p.startswith("Project_")
        IDs.append(p.split("_")[1])
        PI = projectPI(p)
    # Internal vs external PIs are shipped differently (see wd40's
    # fetchFolders): external PIs only get their fastqs fex'ed and never
    # get an internal sequencing_data directory, so there is nothing here
    # for this tool to point users at yet.
    if forcePI is not None:
        PI = forcePI.lower()
    elif PI not in config["Internals"]["PIs"].split(","):
        sys.exit(
            f"PI '{PI}' is not in the internal PI list, so this project was "
            "likely delivered externally via Fex (same check as 'wd40 rel .' "
            "does), which doesn't produce an internal sequencing_data "
            "directory. emailProjectFinished doesn't support externally "
            "fex'ed projects yet."
        )
    # Get the actual sequencing_data dir
    # Assume if multiple projects are given, they all in the same flowcell.
    flowcell = getFlowCell()
    # Assume that only a flow cell exists only once.
    matches = glob.glob(
        os.path.join(
            config["Dirs"]["piDir"], PI, config["Internals"]["seqDir"] + "*", flowcell
        )
    )
    prefix = config["Internals"]["seqDir"]

    def volume_number(match):
        suffix = os.path.basename(os.path.dirname(match))[len(prefix) :]
        if suffix.startswith("_"):
            suffix = suffix[1:]
        if not suffix:
            return 0
        return int(suffix) if suffix.isdecimal() else -1

    matches = [match for match in matches if volume_number(match) >= 0]
    if not matches:
        sys.exit(
            f"No sequencing_data directory found for PI '{PI}' and flowcell "
            f"'{flowcell}' under {config['Dirs']['piDir']}. Double check the "
            "project was actually shipped internally."
        )
    selected = max(matches, key=volume_number)
    seqdir = os.path.basename(os.path.dirname(selected))

    if len(IDs) == 1:
        return IDs[0], seqdir

    return " and ".join([", ".join(IDs[:-1]), IDs[-1]]), seqdir


def getFlowCell():
    return os.path.split(os.getcwd())[-1]


def parse_force_to(value):
    parts = value.split(",")
    if len(parts) != 2:
        raise ValueError("--force-to must be EMAIL,PI")
    email, PI = (part.strip() for part in parts)
    if not email or "@" not in email or any(char in email for char in "\r\n"):
        raise ValueError("--force-to must contain an email address")
    if not PI or any(separator in PI for separator in ("/", "\\", ",")):
        raise ValueError("--force-to must contain a simple PI name")
    return email, PI.lower()


def _run_fexsend(args, allow_failure=False):
    try:
        output = sp.check_output(["fexsend", *args], stderr=sp.STDOUT)
    except sp.CalledProcessError as e:
        if not allow_failure:
            raise RuntimeError(
                f"fexsend {' '.join(args)} failed with exit code {e.returncode}"
            ) from e
        output = e.output or b""
    except OSError as e:
        raise RuntimeError(f"Unable to run fexsend: {e}") from e
    if isinstance(output, bytes):
        return output.decode("utf-8", errors="replace")
    return output


def _fex_archive_names(project):
    return [
        f"{getFlowCell()}_{project}_ro_crate.zip",
        f"{project}_rocrate.zip",
    ]


def _fex_list_entry(output, archive_names):
    entry = None
    for line in output.splitlines():
        for archive_name in archive_names:
            if archive_name in line:
                match = re.search(r"#?\s*(\d+)\)", line)
                entry = archive_name, int(match.group(1)) if match else None
    return entry


def _fex_link_from_output(output, archive_name, allow_any=False):
    for match in re.finditer(r"https?://[^\s<>\"']+", output):
        url = match.group(0).rstrip(".,;)")
        if not allow_any and archive_name not in unquote(url):
            continue
        parsed = urlsplit(url)
        path_parts = [part for part in parsed.path.split("/") if part]
        if "fop" not in path_parts:
            continue
        fop_index = path_parts.index("fop")
        if len(path_parts) <= fop_index + 1:
            continue
        dkey = path_parts[fop_index + 1]
        if not dkey:
            continue
        filename_parts = path_parts[fop_index + 2 :]
        if filename_parts and filename_parts[-1] not in {dkey, "LIST"}:
            path = parsed.path
        else:
            prefix = "/".join(path_parts[:fop_index])
            path = f"{prefix}/fop/{quote(dkey, safe='')}/{quote(archive_name, safe='')}"
        return urlunsplit((parsed.scheme, parsed.netloc, path, "", ""))
    return None


def getFexLinks(projects, config):
    from_address = config["communication"]["fromaddress"]
    list_output = _run_fexsend(["-l", from_address])
    links = {}
    for project in projects:
        project_name = os.path.basename(os.path.normpath(project))
        archive_names = _fex_archive_names(project_name)
        entry = _fex_list_entry(list_output, archive_names)
        if entry is None:
            archive_list = ", ".join(archive_names)
            raise RuntimeError(
                f"Could not find FEX archive for {project_name} using "
                f"fexsend -l; checked {archive_list}"
            )
        archive_name, file_number = entry
        link = _fex_link_from_output(list_output, archive_name)
        if link is None and file_number is not None:
            detail_output = _run_fexsend(
                ["-l", str(file_number), from_address], allow_failure=True
            )
            link = _fex_link_from_output(detail_output, archive_name, allow_any=True)
        if link is None:
            raise RuntimeError(
                f"Could not determine the FEX download link for {archive_name}"
            )
        links[project] = link
    return links


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Send an email to one or more users about a project(s) \
             being finished. This must be run in the output directory of the \
            demultiplexing pipeline."
    )
    parser.add_argument(
        "--configfile",
        default=os.path.expanduser("~/configs/dissectBCL_prod.ini"),
        help="specify a custom ini file. default = {}".format(
            os.path.expanduser("~/configs/dissectBCL_prod.ini")
        ),
    )
    parser.add_argument(
        "--notGood",
        action="store_true",
        help="If specified, \
        do NOT say that the sequencing quality was good.",
    )
    parser.add_argument(
        "--analysis",
        action="store_true",
        help="If specified, \
        the BigRedButton did something with these projects.",
    )
    parser.add_argument("--cc", nargs="+", help="One or more addresses to CC.")
    parser.add_argument(
        "--comment",
        help="Either comment that will be \
        included as its own paragraph (ensure you quote the whole thing!) or \
        the path to a file containing such a comment.",
    )
    parser.add_argument(
        "--fromPerson", help="The name of the person sending the email."
    )
    parser.add_argument(
        "--fromEmail",
        help="The email address of the person \
        sending this. Note that they receive a copy as BCC!",
    )
    parser.add_argument(
        "--fromSignature",
        help="An optional signature of the person \
        sending this.",
    )
    parser.add_argument(
        "--toEmail",
        help="The email address of the person \
         who will receive this.",
        default="",
    )
    parser.add_argument(
        "--force-to",
        metavar="EMAIL,PI",
        help="Force the recipient and sequencing-data PI, and add the FEX download link to the comments.",
    )
    parser.add_argument(
        "--toName",
        help="The name of the person who will \
        receive this.",
        default="",
    )
    parser.add_argument(
        "project",
        nargs="+",
        help="One or more project \
        directories. Only the user on the first will receive an email!",
    )

    args = parser.parse_args(argv)
    force_to = None
    if args.force_to is not None:
        try:
            force_to = parse_force_to(args.force_to)
        except ValueError as e:
            parser.error(str(e))

    print(f"emailProjectFinished: Loading conf from {args.configfile}")
    config = getConf(args.configfile, quickload=True)

    # Double check the project folder(s) actually exist.
    for p in args.project:
        if not os.path.exists(p):
            sys.exit(f"Project folder {p} not found.")

    if force_to:
        try:
            fexLinks = getFexLinks(args.project, config)
        except RuntimeError as e:
            sys.exit(f"emailProjectFinished: {e}")
    else:
        fexLinks = {}

    # get user from project name, lastName = args.project[0].split("_")[2]
    if force_to:
        firstName = args.toName or "there"
        email = force_to[0]
    elif not args.toEmail or not args.toName:
        my_dict = getContactDetails(args.project[0].split("_")[1], config)
        firstName, email = my_dict["first_name"], my_dict["email"]
    else:
        firstName, email = args.toName, args.toEmail

    if not firstName or not email:
        sys.exit("User is not known or does not have an email!\n")

    if not args.fromPerson or not args.fromEmail:
        sys.exit("Sender is not known or does not have an email!\n")

    content = f"""Hi {firstName},

Your sequencing samples for project"""

    if len(args.project) > 1:
        content += "s"
    if force_to:
        project_ids, seqdir = getProjectIDs(args.project, config, forcePI=force_to[1])
    else:
        project_ids, seqdir = getProjectIDs(args.project, config)
    content += (
        f" {project_ids} are finished and the results are now available in your "
        f"group's {seqdir} directory"
    )

    content += f" under the {getFlowCell()} folder.\n"

    if not args.notGood:
        content += "The overall sequencing quality for these samples was good."

    if args.analysis:
        content += (
            "\nAn automated partial analysis "
            "(https://doi.org/10.1093/bioinformatics/btz436) "
            "is available in the same location. \nIf you would like completed "
            "analysis in a semi-automated fashion, please request that using "
            "our online portal: http://snakequest.ie-freiburg.mpg.de .\n"
        )

    comments = []
    if args.comment:
        if os.path.exists(args.comment):
            with open(args.comment) as commentFile:
                comments.append(commentFile.read())
        else:
            comments.append(args.comment)
    for project, link in fexLinks.items():
        comments.append(f"FEX download link for {project}: {link}")
    if comments:
        content += "\n===\n" + "\n\n".join(comments) + "\n===\n"

    content += f"\nPlease let me know if you have any other questions,\
        \n{args.fromPerson}\n"

    # Add a .signature
    if args.fromSignature is not None and os.path.exists(args.fromSignature):
        content += "\n--\n"
        with open(args.fromSignature) as signatureFile:
            content += signatureFile.read()

    # Send the Email
    msg = MIMEText(content)
    msg["Subject"] = "Sequencing samples ready - " + args.project[0]
    msg["From"] = args.fromEmail
    msg["To"] = email
    bioinfo_cc = config["communication"]["bioinfoCore"]
    host = config["communication"]["host"]

    if args.cc:
        msg["Cc"] = ", ".join(args.cc)

    msg["Bcc"] = bioinfo_cc

    s = smtplib.SMTP(host)
    s.send_message(msg)
    s.quit()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("--help")
    main()
