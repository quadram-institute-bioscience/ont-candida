#!/usr/bin/env python
"""
Andrea Telatin 2021
Align two sequences with lastz and plot the results with gnuplot.

REQUIRES:
    lastz
    gnuplot / R
"""

import copy
import sys, os
import subprocess
import argparse
import tempfile
from dataclasses import dataclass # Python 3.7+
import shutil
import logging

# Make logger
logging.basicConfig(
        format='%(levelname)-8s :%(asctime)s: %(message)s', 
        datefmt='%H:%M:%S') #%(name)s
logger = logging.getLogger(__name__)

@dataclass
class LastzOpts:
    format: str = "rdotplot"
    strand: str = "both"
    step: int = 20
    nogapped: bool = True
    notransition: bool = True


    def asString(self):
        return " ".join(self.asList())

    def asList(self):
        opts = []
        opts.append("--format={}".format(self.format))
        opts.append("--strand={}".format(self.strand))
        opts.append("--step={}".format(self.step))
        if self.nogapped:
            opts.append("--nogapped")
        if self.notransition:
            opts.append("--notransition")
        return opts


def check_dependency(cmd):
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        logger.info("Found dependency: {}".format(cmd[0]))
        return True
    except subprocess.CalledProcessError as e:
        logger.error("Missing dependency: {}".format(cmd[0]))
        return False

def fileExistsNotEmpty(filename):
    return os.path.isfile(filename) and os.stat(filename).st_size > 0

def lastz(referenceSeq, querySeq, outfile, options):
    """
    Run LastZ on the two sequences
    """
    cmd = ["lastz"] + options.asList() + [referenceSeq, querySeq]
    logger.info("Running LastZ: {}".format(outfile))
    logger.debug("CMD: " + " ".join(cmd))
    with open(outfile, "w") as f:
        try:
            run = subprocess.run(cmd, stdout=f)
            if not fileExistsNotEmpty(outfile):
                raise Exception("lastz failed: {} is empty".format(outfile))
            return run.returncode
        except subprocess.CalledProcessError as e:
            print("Error executing LastZ: {}".format(e))
            sys.exit(1)
        except Exception as e:
            print("Unexpected LastZ error: {}".format(e))
            sys.exit(1)


def rdotplot(dotfile, output, format="pdf"):
    script = """
    plotData <- read.delim("{dotfile}");
    {format}("{output}");
    plot(plotData, type="l");
    dev.off();
    """.format(dotfile=os.path.abspath(dotfile), output=os.path.abspath(output), format=format)
    cmd = [ "Rscript", "--vanilla", "--no-save", "-e", script]
    try:
        logger.info("running: {}".format(" ".join(cmd)))
        run = subprocess.run(cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        if not fileExistsNotEmpty(output):
            raise Exception("Rscript failed: {} is empty".format(output))
        return run.returncode
    except subprocess.CalledProcessError as e:
        logger.error("Error executing Rscript: {}".format(e))
        sys.exit(1)

def gnuplot_pipeline(seq1, seq2, lastz_output, out_gnuplot, out_gnuplotscript, out_png):
    try:
        logger.info("[GNUPLOT] Converting rplot to gnuplot: {}".format(out_gnuplot))
        gnuplot_r2gnu(lastz_output, out_gnuplot)
    except Exception as e:
        print("Error running `gnuplot_r2gnu`: {}".format(e), file=sys.stderr)
        sys.exit(1)

    try:
        logger.info("[GNUPLOT] Generating GNU plot script: {}".format(out_gnuplotscript))
        gnuplot_makescript(seq1, seq2, out_gnuplot, out_gnuplotscript, out_png)
    except Exception as e:
        print("Error running `gnuplot_makescript`: {}".format(e), file=sys.stderr)
        sys.exit(1)

    try:
        logger.info("[GNUPLOT] Generating GNU plot: {}".format(out_png))
        gnuplot_runscript(out_gnuplotscript)

    except Exception as e:
        print("Error running `gnuplot_runscript`: {}".format(e), file=sys.stderr)
        sys.exit(1)
def gnuplot_r2gnu(input, gnuplot):
    """
    Convert a rplot file to a gnuplot-readable file.
    """
    #cat seq1-23.rdotplot|  paste - - - | awk '{printf("%d %d %d %d\n",$3,$4,int($5)-int($3),int($6)-int($4));}' > seq1-23.gnuplot
    with open(input, "r") as f:
        with open(gnuplot, "w") as g:
            c = 0
            fields = []
            sum1 = 0
            sum2 = 0
            for line in f:
                c += 1
                f = line.strip().split("\t")
                if c == 1:
                    seq1 = f[0]
                    seq2 = f[1]
                fields = fields + f

                # every three lines
                if c % 3 == 0:
                    #print(fields)
                    try:
                        output = [fields[2], fields[3],int(fields[4]) - int(fields[2]), int(fields[5]) - int(fields[3])]
                        sum1 += output[2]
                        sum2 += output[3]
                        g.writelines("\t".join(map(str, output)) + "\n")

                    except Exception as e:
                        print("Error: {}".format(e))
                        print("2:", fields[2])
                        sys.exit(1)

                    fields.clear()
                
                #g.write("{} {} {} {}\n".format(fields[3], fields[4], int(fields[5])-int(fields[3]), int(fields[6])-int(fields[4])))
    #print("Lengths: {}, {}".format(sum1, sum2))

def gnuplot_makescript(seq1, seq2, datafile, script, plotfile):
    seq2name = os.path.basename(seq2).split(".")[0].replace("_", " ").replace("-", " ")
    seq1name = os.path.basename(seq1).split(".")[0].replace("_", " ").replace("-", " ")
    gnuScriptText = """set terminal png;
set xlabel '{seq2}';
set ylabel '{seq1}';
set output "{outfile}";
set xtics rotate;
plot "{script}" using 1:2:3:4 with vectors nohead
            """.format(seq1=seq1name, seq2=seq2name, outfile=plotfile, script=os.path.basename(datafile))
    with open(script, "w") as f:
        f.write(gnuScriptText)

def gnuplot_runscript(script):
    # Change directory to the script's directory
    currentdir = os.getcwd()
    os.chdir(os.path.dirname(script))
    cmd = ["gnuplot", script]
    try:
        run = subprocess.run(cmd)
        os.chdir(currentdir)
        return run.returncode
    except subprocess.CalledProcessError as e:
        print("Error executing gnuplot: {}".format(e))
        os.chdir(currentdir)
        sys.exit(1)

def parseAln(alnfile):
    ##nmismatch      name1   strand1 start1  end1    name2   strand2 start2  end2
    #2       contig_22       +       202     405     Contig005504_C_parapsilosis_CDC317      +       1       204
    mismatches = 0
    query = 0
    target = 0

    with open(alnfile, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            mismatches += int(fields[0])
            target += int(fields[4]) - int(fields[3])
            query += int(fields[8]) - int(fields[7])
    return mismatches, query, target
if __name__ == "__main__":
    defaultTempDir = tempfile.gettempdir()
    args = argparse.ArgumentParser(description="Align two sequences with lastz and plot the results with gnuplot.")
    args.add_argument("DB", help="Reference database in FASTA format")
    args.add_argument("QUERY", help="Query sequence in FASTA format")
    args.add_argument("-o", "--basename", help="Output file basename", required=True)
    args.add_argument("-n", "--plotname", help="Plot file name [default: %(default)s", default="Contig_alignemnt")
    args.add_argument("--tmp", help="Where to write temp directory (must exist) [default: %(default)s", default=defaultTempDir)
    args.add_argument("--keeptmp", help="Keep temporary directory", action="store_true")
    args.add_argument("--gnuplot", help="Use GNU plot instead of R", action="store_true")
    args.add_argument("--pdf", help="Save in PDF format (only with R)", action="store_true")
    
    args.add_argument("--verbose", help="Print verbose output", action="store_true")
    args = args.parse_args()

    if args.verbose:
        LEVEL = logging.DEBUG
    else:
        LEVEL = logging.WARNING



    
    logger.setLevel(LEVEL)

    # Create temporary files
    outformat = "png" if args.gnuplot else "pdf" if args.pdf else "png"
    logger.info("Output extension: {}".format(outformat))
    tmp_dir = tempfile.mkdtemp(prefix="lastzPlot_", dir=args.tmp)
    tmp_dir = os.path.abspath(tmp_dir)
    lastz_tmpfile = os.path.join(tmp_dir, "lastz.rdotplot")
    lastz_tmpaln = os.path.join(tmp_dir, "lastz.aln")
    gnuplot_data = os.path.join(tmp_dir, args.plotname)
    gnuscript_tmp = os.path.join(tmp_dir, "Script.gnu")
    plotfile_tmp = os.path.join(tmp_dir, "Plot_Alignment.{ext}".format(ext=outformat))
    lastZoptions = LastzOpts(
        format="rdotplot",
        strand="both",
        step=20,
        nogapped=True,
        notransition=True
    )
    
    logger.info("Temporary directory: {}".format(tmp_dir))
    logger.info("Temporary file: {}".format(os.path.basename(lastz_tmpfile)))
    logger.info("Temporary gnuplot file: {}".format(os.path.basename(gnuplot_data)))
    logger.info("LastZ options: {}".format(lastZoptions.asList()))

    
    # Align
    lastz(args.DB, args.QUERY, lastz_tmpfile, lastZoptions)
    alnOpt = copy.deepcopy(lastZoptions)
    alnOpt.format = "general:nmismatch,name1,strand1,start1,end1,name2,strand2,start2,end2"
    logger.debug(alnOpt.asList())
    
    lastz(args.QUERY, args.DB, lastz_tmpaln, alnOpt)
    mm, q, t = parseAln(lastz_tmpaln)
    logger.info("Aligned: {:,.0f} query, {:,.0f} target [{:,.0f} mismatches]".format(q, t, mm))
    

    if args.gnuplot:
        # PLOT WITH GNUPLOT
        check_dependency(["gnuplot", "--version"])
        gnuplot_pipeline(args.DB, args.QUERY, lastz_tmpfile, gnuplot_data, gnuscript_tmp, plotfile_tmp)
    else:
        # PLOT WITH R
        check_dependency(["Rscript", "--version"])
        rdotplot(lastz_tmpfile, plotfile_tmp, outformat)
  
    # Move plot to the final output file
    if os.path.exists(plotfile_tmp):
        logger.info("Moving {} to {}".format(plotfile_tmp, args.basename + ".{ext}".format(ext=outformat)))
        os.rename(plotfile_tmp, args.basename + ".{ext}".format(ext=outformat))
        logger.info("Plotting result: {}".format(args.basename + ".{ext}".format(ext=outformat)))
    else:
        logger.error("Plotting failed: {} not found".format(plotfile_tmp))

    # Remove temporary files if needed
    if not args.keeptmp:
        try:
            shutil.rmtree(tmp_dir)
            logger.info("Removing temporary directory: {}".format(tmp_dir))
        except Exception as e:
            logger.error("Error removing temporary directory: {}".format(e))
            exit(1)

    