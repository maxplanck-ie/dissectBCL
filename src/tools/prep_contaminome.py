import rich_click as click
from rich import print
import os
import yaml
import shutil
import gzip
import urllib.request
import subprocess
import glob

ignore_chrs = {
    'human': ['NC_012920.1'],  # human mito
    'mouse': ['NC_005089.1'],  # mouse mito
    'fly': ['NC_024511.2'],  # fly mito
    'aedes-aegypti': ['NC_035159.1'],  # aedes mito
    'zebrafish': ['NC_002333.2']  # zebrafish mito
}

rrna_mask = [
    ('human', 'humanrrna'),
    ('mouse', 'mouserrna'),
    ('fly', 'flyrrna'),
    ('aedes-aegypti', 'aedesaegyptirrna'),
    ('zebrafish', 'zebrafishrrna')
]

taxmap = {
    'root': [1, 1, 'no rank'],
    'alive': [2, 1, 'no rank'],
    'non-alive': [3, 1, 'no rank'],
    'eukaryote': [4, 2, 'domain'],
    'prokaryote': [5, 2, 'domain'],
    'phage': [6, 3, 'domain'],
    'virus': [7, 3, 'domain'],
    'vector': [8, 3, 'domain'],
    'humangrp': [9, 4, 'family'],
    'mousegrp': [10, 4, 'family'],
    'flygrp': [11, 4, 'family'],
    'pseudomonasgrp': [12, 5, 'family'],
    'eugrp': [13, 4, 'family'],
    'progrp': [14, 5, 'family'],
    'human': [9606, 9, 'species'],
    'mouse': [10090, 10, 'species'],
    'fly': [7227, 11, 'species'],
    'aedes-aegypti': [7159, 13, 'species'],
    'sea-lamprey': [7757, 13, 'species'],
    'japanese-medaka': [8090, 13, 'species'],
    'c-elegans': [6239, 13, 'species'],
    'yeast': [4932, 13, 'species'],
    'fission-yeast': [4896, 13, 'species'],
    'zebrafish': [7955, 13, 'species'],
    'moss-piglet': [232323, 13, 'species'],
    'ecoli': [562, 14, 'species'],
    'pseudomonas-aeruginosa': [287, 12, 'species'],
    'pseudomonas-putidia': [1211579, 12, 'species'],
    'pseudomonas-fulva': [47880, 12, 'species'],
    'pseudomonas-fluorescens': [294, 12, 'species'],
    'pseudomonas-yamanorum': [515393, 12, 'species'],
    'mycoplasma': [2148, 14, 'species'],
    'mycoplasma-hyorhinis': [2100, 14, 'species'],
    'haemophilus': [727, 14, 'species'],
    'wolbachia': [163164, 14, 'species'],
    'staph_aureus': [1280, 14, 'species'],
    'burkholderia-ubonensis': [1249668, 14, 'species'],
    'sars-cov2': [2697049, 7, 'species'],
    'flu-a': [335341, 7, 'species'],
    'noro': [11983, 7, 'species'],
    'flu-b': [11520, 7, 'species'],
    'common-cold-a': [573824, 7, 'species'],
    'common-cold-b': [12131, 7, 'species'],
    'drosophila-c-virus': [64279, 7, 'species'],
    'phix': [2886930, 6, 'species'],
    'lambdaphage': [10710, 6, 'species'],
    'vectors': [29278, 8, 'species'],
    'ercc': [292781111, 8, 'species'],
    'humanrrna': [96061111, 9, 'species'],
    'mouserrna': [100901111, 10, 'species'],
    'aedesaegyptirrna': [71591111, 13, 'species'],
    'zebrafishrrna': [79551111, 13, 'species'],
    'flyrrna': [72271111, 11, 'species'],
    'humanmito': [96062222, 9, 'species'],
    'mousemito': [100902222, 10, 'species'],
    'flymito': [72272222, 11, 'species'],
    'aedesaegyptimito': [71592222, 13, 'species'],
    'zebrafishmito': [79552222, 13, 'species']
    }


def dumpnode(taxmap, ofile):
    with open(ofile, 'w') as f:
        for taxon in taxmap:
            f.write(
                '{}\t|\t{}\t|\t{}\t|\t-\t|\n'.format(
                    taxmap[taxon][0],
                    taxmap[taxon][1],
                    taxmap[taxon][2],
                )
            )


def dumpname(taxmap, ofile):
    with open(ofile, 'w') as f:
        for taxon in taxmap:
            f.write(
                '{}\t|\t{}\t|\t-\t|\tscientific name\t|\n'.format(
                    taxmap[taxon][0],
                    taxon
                )
            )


def mask(mtupe, libdir):
    '''
    mask sequence mtupe[1] in mptupe[0]
    '''
    print("Building blast db {}".format(mtupe))
    bldb = [
        'makeblastdb',
        '-dbtype',
        'nucl',
        '-in',
        os.path.join(libdir, mtupe[0] + '.fna'),
        '-out',
        'tmpdb'
    ]
    print(bldb)
    subprocess.run(bldb)
    bln = [
        'blastn',
        '-query',
        os.path.join(libdir, mtupe[1] + '.fna'),
        '-out',
        'tmpdb_bln',
        '-db',
        'tmpdb',
        '-evalue',
        '1e-10',
        '-outfmt',
        '6'
    ]
    print(bln)
    subprocess.run(bln)
    bed = []
    with open('tmpdb_bln') as f:
        for line in f:
            lis = line.strip().split()
            if int(lis[8]) > int(lis[9]):
                bed.append(
                    [lis[1], lis[9], lis[8]]
                )
            else:
                bed.append(
                    [lis[1], lis[8], lis[9]]
                )
    with open('tmpdb_bed', 'w') as f:
        for b in bed:
            f.write(
                "{}\t{}\t{}\n".format(b[0], b[1], b[2])
            )
    # mask fna
    fnain = os.path.join(libdir, mtupe[0] + '.fna')
    fnaout = os.path.join(libdir, mtupe[0] + 'masked.fna')
    maskcmd = [
        'bedtools',
        'maskfasta',
        '-fi',
        fnain,
        '-bed',
        'tmpdb_bed',
        '-fo',
        fnaout
    ]
    print(maskcmd)
    subprocess.run(maskcmd)
    # Remove original fasta.
    os.remove(fnain)
    os.rename(fnaout, fnain)
    for tmpf in glob.glob('tmpdb*'):
        os.remove(tmpf)


@click.command(
    context_settings=dict(
        help_option_names=["-h", "--help"]
    )
)
@click.option(
    '-c',
    '--contaminome',
    required=True,
    type=click.Path(exists=True),
    help='Specify contaminome.yml file'
)
@click.option(
   "-o",
   "--outputdir",
   type=click.Path(exists=True),
   required=True,
   help="Specify an output directory."
)
@click.option(
   "-t",
   "--threads",
   required=False,
   default=15,
   help="Set number of threads"
)
def main(contaminome, outputdir, threads):
    contaminomedir = os.path.join(
        outputdir,
        'contaminomedb'
    )
    # Purge existing contaminome folder.
    if os.path.exists(contaminomedir):
        print("{} already exists. Purging and re-creating.".format(
            contaminomedir
        ))
        shutil.rmtree(
            contaminomedir
        )
    os.mkdir(contaminomedir)
    librarydir = os.path.join(
        contaminomedir,
        'library'
    )
    os.mkdir(librarydir)
    # read up yaml file.
    with open(contaminome) as f:
        genomes = yaml.safe_load(f)

    seq2acc = {}
    # Download, format
    for domain in genomes:
        for sp in genomes[domain]:
            print("Fetching - {}".format(sp))
            if genomes[domain][sp]['vulgarname'] in list(
                ignore_chrs.keys()
            ):
                ignoring = ignore_chrs[
                    genomes[domain][sp]['vulgarname']
                ]
            else:
                ignoring = []
            ofile = os.path.join(
                librarydir,
                genomes[domain][sp]['vulgarname'] + '.fna'
            )
            if not os.path.exists(ofile):
                with urllib.request.urlopen(
                    genomes[domain][sp]['URL']
                ) as gzip_genome:
                    with gzip.GzipFile(fileobj=gzip_genome) as f:
                        genome = f.read().decode().splitlines()
                headCount = 0
                with open(ofile, 'w') as f:
                    for line in genome:
                        if line.startswith('>'):
                            headCount += 1
                            headStr = line.strip().replace(
                                '>', ''
                            ).replace(
                                '|', ''
                            ).replace(
                                'gnluv',
                                ''
                            ).replace(
                                ':',
                                ''
                            ).replace(
                                '-',
                                ''
                            ).split(' ')[0]
                            # Check we need to ignore this chromosome
                            if headStr in ignoring:
                                appendStatus = False
                            else:
                                appendStatus = True
                            if appendStatus:
                                headStr += '|kraken:taxid|'
                                headStr += str(
                                    genomes[domain][sp]['taxid']
                                )
                                f.write('>' + headStr + '\n')
                                seq2acc[headStr] = genomes[
                                    domain
                                ][sp]['taxid']
                        else:
                            if appendStatus:
                                f.write(line.strip() + '\n')
    seq2taxfile = os.path.join(contaminomedir, 'seqid2taxid.map')
    with open(seq2taxfile, 'w') as f:
        for i in seq2acc:
            f.write("{}\t{}\n".format(
                i, seq2acc[i]
            ))
    taxondir = os.path.join(contaminomedir, 'taxonomy')
    os.mkdir(taxondir)
    nodesfile = os.path.join(taxondir, 'nodes.dmp')
    namesfile = os.path.join(taxondir, 'names.dmp')
    dumpnode(taxmap, nodesfile)
    dumpname(taxmap, namesfile)
    # Mask rRNA if we need to.
    for masktupe in rrna_mask:
        mask(masktupe, librarydir)

    # Build kraken db
    krakcmd = [
        'kraken2-build',
        '--build',
        '--threads',
        str(threads),
        '--db',
        contaminomedir
    ]
    subprocess.run(krakcmd)
