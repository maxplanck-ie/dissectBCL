#!/package/anaconda3/bin/python3
import argparse
import sys
import smtplib
import os
import smtplib
from email.mime.text import MIMEText


def fetchFirstNameAndEmail(lastName):
    # search in dictionary file for lastName
    fn="/home/pipegrp/configs/parkourUsers.txt"
    f = open(fn)
    d=dict() 
    for line in f:
        cols = line.rstrip().split("\t")

        # only accept format: firstName, lastName, email
        if (len(cols) < 3):
            continue
   
        # ignore all other lastNames
        if cols[1] != lastName:
            continue

        # check if lastName occurs more than once in list
        if cols[1] in d:
            print("Error: fetchFirstNameAndEmail\n\
Name {} exists more than once. Specify --toEmail and --toName explicitly!".format(cols[1]))
            print('now:      ', cols[1], cols[0], cols[2])
            print('previous: ', cols[1], d[cols[1]])
            sys.exit(1)

        # add to dictionary
        d[cols[1]] = [cols[0], cols[2]]
    f.close()

    if lastName not in d:
        print("Error: fetchFirstNameAndEmail\n\
No Information for lastName={}. You may have to update {}".format(lastName, fn))
        sys.exit(1)

    return d[lastName]


def getProjectIDs(projects):
    IDs = []
    for p in projects:
        # Sanity check
        assert(p.startswith("Project_"))
        IDs.append(p.split("_")[1])
    if len(IDs) == 1:
        return IDs[0]
    return " and ".join([", ".join(IDs[:-1]), IDs[-1]])


def getFlowCell():
    return os.path.split(os.getcwd())[-1]


def main():
    parser = argparse.ArgumentParser(description="Send an email to one or more users about a project or project being finished. This must be run in the output directory of the demultiplexing pipeline.")
    parser.add_argument("--notGood", action="store_true", help="If specified, do NOT say that the sequencing quality was good.")
    parser.add_argument("--analysis", action="store_true", help="If specified, the BigRedButton did something with these projects.")
    parser.add_argument("--cc", nargs="+", help="One or more addresses to CC.")
    parser.add_argument("--comment", help="Either comment that will be included as its own paragraph (ensure you quote the whole thing!) or the path to a file containing such a comment.")
    parser.add_argument("--noGalaxy", action="store_true", help="Do NOT say that files are in Galaxy.")
    parser.add_argument("--fromPerson", help="The name of the person sending the email. Note that the contents of ~/.signature_fromPerson are also used! (Default: manke)", default="manke")
    parser.add_argument("--fromEmail", help="The email address of the person sending this. Note that they receive a copy as BCC!", default="manke@ie-freiburg.mpg.de")
    parser.add_argument("--toEmail", help="The email address of the person who will receive this.", default="")
    parser.add_argument("--toName",  help="The name of the person who will receive this.", default="")
    parser.add_argument("project", nargs="+", help="One or more project directories. Only the user on the first will receive an email!")
    args = parser.parse_args()

    # get lastName (user) from project name
    lastName = args.project[0].split("_")[2]
    if not args.toEmail or not args.toName:
        firstName, email = fetchFirstNameAndEmail(lastName)
    else:
        firstName, email = args.toName, args.toEmail

    if not firstName or not email:
        sys.exit("User is not known or does not have an email!\n")

    content = """Hi {},

Your sequencing samples for project""".format(firstName)

    if len(args.project) > 1:
        content += "s"
    content += " {} are finished and the results are now available in your group's sequencing data directory".format(getProjectIDs(args.project))

#    Stop saying galaxy alltogether, We stop automatically uploading.
#    if not args.noGalaxy:
#        content += " and Galaxy"
    content += " under the {} folder.".format(getFlowCell())

    if not args.notGood:
        content += " The overall sequencing quality for these samples was good."
    if args.analysis:
        content += " An automated partial analysis (https://doi.org/10.1093/bioinformatics/btz436) of these samples is available in the same location. If you would like completed analysis of these samples in a semi-automated fashion, please request that using our online portal: http://snakequest.ie-freiburg.mpg.de .\n"

    content += " Please note that sequencing data is no longer deposited into Galaxy by default. If you need to access this data in Galaxy, please let me know. \n"

    if args.comment:
        content += "\n"
        if os.path.exists(args.comment):
            content += open(args.comment).read()
        else:
            content += args.comment
        content += "\n"

    content += "\nPlease let me know if you have any other questions,\n{}\n".format(args.fromPerson)

    # Add a .signature
    if os.path.exists("/home/pipegrp/.signature_"+args.fromPerson):
        content += "\n--\n"
        content += open("/home/pipegrp/.signature_"+args.fromPerson).read()

    # Send the Email
    msg = MIMEText(content)
    msg['Subject'] = "Sequencing samples ready - " + args.project[0]
    msg['From'] = args.fromEmail
    msg['To'] = email
    if args.cc:
#        args.cc.append("bioinfo-core@ie-freiburg.mpg.de")
        msg['Cc'] = ", ".join(args.cc)
    else:
        msg['Cc'] = ''
#        msg['Cc'] = "bioinfo-core@ie-freiburg.mpg.de"
    msg['Bcc'] = args.fromEmail

    s = smtplib.SMTP("mail.ie-freiburg.mpg.de")
    s.send_message(msg)
    s.quit()

if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("--help")
    main(sys.argv[1:])
