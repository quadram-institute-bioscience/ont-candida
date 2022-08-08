#!/usr/bin/env python
import copy
from os import remove
from weakref import ref
from gnuplotter import *
import re
def read_fasta(path):
    if path.endswith(".gz"):
        import gzip
    name = None
    with (gzip.open if path.endswith('.gz') else open)(path, 'rt') as fasta:
        for line in fasta:
            if line.startswith('>'):
                if name is not None:
                    yield name, seq
                name = line[1:].rstrip().split()[0]
                seq = ''
            else:
                seq += line.rstrip()
    yield name, seq

def stringToFilename(s):
    # Remove all non-alphanumeric characters
    s = re.sub('[^a-zA-Z0-9]', '', s)
    # Remove all leading and trailing whitespace
    s = s.strip()
    # Replace whitespace with underscores
    s = s.replace(' ', '_')
    return s

def split_fasta(file, destdir, tag):
    files = {}
    c = 0
    
    for name, seq in read_fasta(file):
        c += 1
        cFmt = str(c).zfill(9)
        filename = '%s_%s.fasta' % (tag, cFmt)
        output = os.path.abspath(os.path.join(destdir, filename))
        files[output] = (name, len(seq))

        with open(output, 'wt') as f:
            f.write('>' + name + '\n')
            f.write(seq)
            f.close()
    return files
if __name__=="__main__":
    args = argparse.ArgumentParser()
    args.add_argument("REF", help="Reference file")
    args.add_argument("QUERY", help="Query file")
    args.add_argument("-o", "--output", help="Output directory", required=True)
    args.add_argument("--verbose", help="Verbose output", action="store_true")
    args.add_argument("--debug", help="Debug output", action="store_true")

    args = args.parse_args()

    if args.debug == True:
        logger.setLevel(logging.DEBUG)
    elif args.verbose == True:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    # Check input files exist
    if not os.path.isfile(args.REF):
        logger.error("Reference file does not exist: %s" % args.REF)
        sys.exit(1)
    
    if not os.path.isfile(args.QUERY):
        logger.error("Query file does not exist: %s" % args.QUERY)
        sys.exit(1)

    tempdir = tempfile.mkdtemp(prefix="compareAndPlot_")
    logger.warning("Tempdir: %s" % tempdir)

    query_file = os.path.join(tempdir, "query.fa")
    ref_file = os.path.join(tempdir, "ref.fa")

    query_files = split_fasta(args.QUERY, tempdir, "query")
    ref_files = split_fasta(args.REF, tempdir, "ref")

    lastzOpts =  LastzOpts(
        format="general:nmismatch,name1,strand1,start1,end1,name2,strand2,start2,end2",
        strand="both",
        step=20,
        nogapped=True,
        notransition=True
    )

    plotOpts = copy.deepcopy(lastzOpts)
    plotOpts.format = "rdotplot"

    refMatches = {}
    i = 0
    for refFile, refTuple  in ref_files.items():
        
        refName, refLen = refTuple
        refMatches[refName] = []
        logger.info("Ref file: %s" % refName)
        for queryFile, queryTuple in query_files.items():
            queryName, queryLen = queryTuple
            

            logger.info("  query file: %s" % queryName)
            tmpout = os.path.join(tempdir, "{}_{}.out".format(stringToFilename(refName), stringToFilename(queryName)))
            logger.debug("Lastz to: %s" % tmpout)
            lastz(refFile, queryFile, tmpout, lastzOpts)
            
            mm, q, t = parseAln(tmpout)
            if t == 0:
                logger.debug("No alignments found for %s and %s" % (refName, queryName))
                continue
            ratioT = 100 * float(t) / float(refLen)
            
            ratioQ = 100 * float(q) / float(queryLen)
            if  ratioT > 12:
  
                refMatches[refName].append(queryName)
          
                logger.info("Ratio: %f : %d - %d" % (ratioT, refLen, t) + "| Ratio: %f : %d - %d" % (ratioQ, queryLen, q))
                logger.info("Aligned: {:,.0f} query, {:,.0f} target [{:,.0f} mismatches]".format(q, t, mm))

                logger.debug("Printing PDF to: {}".format(tmpout + ".pdf"))
                lastz(refFile, queryFile, tmpout + ".plot", plotOpts)
                rdotplot(tmpout + ".plot", tmpout + ".pdf", format="pdf")

        print(refName, ": len=", refLen, " matches=", len(refMatches[refName]), " data=", ",".join(refMatches[refName]), sep="")
        for m in refMatches[refName]:
            l = "?"
            for T, TUPLE in query_files.items():
                if TUPLE[0] == m:
                    l = TUPLE[1]
                    break
            i += 1
            print(i, "  ctg=", m, " len=", l, sep="")

