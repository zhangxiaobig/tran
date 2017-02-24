# python modules
import sys
import os
import logging


def read_input(parser):
    args = parser.parse_args()
    if not os.path.isfile(args.tefile) :
        logging.error("No such file: %s !\n" %(args.tefile))
        sys.exit(1)
    if not os.path.isfile(args.gtffile) :
        logging.error("No such file: %s !\n" % (args.gtffile))
        sys.exit(1)
    # Obtain & store list of files for group 1 (e.g. treatment/mutant)
    for i in range(len(args.tfiles)) :
        if not os.path.isfile(args.tfiles[i]) :
            logging.error("No such file: %s !\n" % (args.tfiles[i]))
            sys.exit(1)
    # Obtain & store list of files for group2 (e.g. control/wildtype)
    for i in range(len(args.cfiles)) :
        if not os.path.isfile(args.cfiles[i]) :
            logging.error("No such file: %s !\n" % (args.cfiles[i]))
            sys.exit(1)
    # Identify file format for subsequent processing (parsing)
    if args.format == "BAM" :
        args.parser = "BAM"
    elif args.format == "SAM" :
        args.parser = "SAM"
    else :
        logging.error("Does not support such file format: %s !\n" % (args.format))
        sys.exit(1)
    # What type of RNA-Seq experiment (stranded or not)
    if args.stranded not in ['yes', 'no', 'reverse'] :
        logging.error("Does not support such stranded value: %s !\n" % (args.stranded))
        sys.exit(1)
        
    logging.basicConfig(level=logging.DEBUG,
        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
        datefmt='%a, %d %b %Y %H:%M:%S',
    stream=sys.stderr, filemode="w")
    
    args.error = logging.critical        # function alias
    args.warn = logging.warning
    args.debug = logging.debug
    args.info = logging.info
    
    args.argtxt = "\n".join(("# ARGUMENTS LIST:", \
                "# name = %s" % (args.prj_name), \
                "# treatment files = %s" % (args.tfiles), \
                "# control files = %s" % (args.cfiles), \
                "# GTF file = %s " % (args.gtffile), \
                "# TE file = %s " % (args.tefile), \
                "# stranded = %s " % (args.stranded), \
    ))
    return args 