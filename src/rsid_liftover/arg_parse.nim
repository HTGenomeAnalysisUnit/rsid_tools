import argparse
import strformat
from ./utils import log

var p = newParser("rsID_liftover"):
    help("")
    arg("indf", nargs = -1, help="input table/list containing rsID to liftover")
    option("-o", "--out", help="Output folder or prefix. If prefix the output file(s) names is: {out}_{infileprefix}.{target}.tsv. If not provided output to stdout")
    flag("-e", "--header", help="Input file has header")
    option("-s", "--sep", help="Delimiter for input file", default = some("\t"))
    option("-m", "--mode", help="Mode of running", default = some("liftover"), choices = @["liftover", "rsid"])
    option("-r", "--rsid_column", help="0-based index of column containing rsIDs in the input", default = some("0"))
    option("-x", "--chrom_column", help="0-based index of column containing chromosome in the input", default = some("1"))
    option("-c", "--chrom", help="Only read data for specific chromosomes. Comma-sep list accepted or -1 for all chromosomes", default = some("-1"))
    option("-t", "--target", help="Target genome build for liftover", required = true, choices = @["GRCh37", "GRCh38"])
    flag("-n", "--no_missing", help="Do not include in output rsIDs that can't be liftovered to the target build")
    option("-d", "--map_dir", help="Directory containing map files for liftover", required = true)
    option("-v", "--version", help="dbSNP version", default = some("151"))       

proc parseCmdLine*(args: seq[string]): ref =
    try:
        result = p.parse(args) 
    except ShortCircuit as e:
        if e.flag == "argparse_help":
            echo p.help
            quit QuitSuccess
    except UsageError:
        stderr.writeLine getCurrentExceptionMsg() 
        echo "Use --help for usage information"
        quit QuitSuccess

proc logArgs*(opts: ref) {.discardable.} =
    log("ARG", fmt"Input file: {opts.indf}")
    let out_stream = (if opts.out != "": opts.out else: "stdout")
    log("ARG", fmt"Output: {out_stream}")
    log("ARG", fmt"Delimiter: {opts.sep}")
    log("ARG", fmt"rsID col index: {opts.rsid_column}")
    log("ARG", fmt"Chrom col index: {opts.chrom_column}")
    log("ARG", fmt"Chromosomes: {opts.chrom}")
    log("ARG", fmt"Target build: {opts.target}")
    log("ARG", fmt"Map folder: {opts.map_dir}")
    log("ARG", fmt"dbSNP version: {opts.version}")
    log("ARG", fmt"Mode: {opts.mode}")
