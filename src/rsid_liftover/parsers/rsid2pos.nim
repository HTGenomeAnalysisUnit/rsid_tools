import argparse
import strformat
from ../utils import log
import ../constants

var p = newParser("liftover"):
    help("a set of small utilities to work with rsIDs")
    arg("intables", nargs = -1, help="input table/list containing rsID to liftover")
    option("-o", "--out", help="Output folder or prefix. Output structure is: rsid2pos_{target_build}-{infileprefix}.tsv. If not provided output to stdout")
    flag("-e", "--header", help="Input file has header")
    option("-s", "--sep", help="Delimiter for input file", default = some("\t"))
    option("-i", "--rsid_column", help="0-based index of column containing rsIDs in the input", default = some($RSID_COLIDX))
    option("-x", "--chrom_column", help="0-based index of column containing chromosome in the input. Use -1 if not present", default = some($CHROM_COLIDX))
    option("-c", "--chrom", help="Only process data for specific chromosomes. Comma-sep list accepted or -1 for all chromosomes", default = some("-1"))
    option("-t", "--target", help="Target genome build for liftover", required = true, choices = @["GRCh37", "GRCh38"])
    flag("-n", "--no_missing", help="Do not include in output rsIDs that can't be liftovered to the target build")
    option("-d", "--map_dir", help="Directory containing binaries files generated by make_bin", required = true)
    run:
        log("ARG", fmt"Output: {opts.out}")
        log("ARG", fmt"Input file(s): {opts.intables.len}")
        log("ARG", fmt"Delimiter: {opts.sep}")
        log("ARG", fmt"Header: {opts.header}")
        log("ARG", fmt"rsID col index: {opts.rsid_column}")
        log("ARG", fmt"Chrom col index: {opts.chrom_column}")
        log("ARG", fmt"Chromosomes: {opts.chrom}")
        log("ARG", fmt"Target build: {opts.target}")
        log("ARG", fmt"Map folder: {opts.map_dir}")

proc parseCmdLine*(args: seq[string]): ref =
    try:
        p.run(args)
        result = p.parse(args) 
    except ShortCircuit as e:
        if e.flag == "argparse_help":
            echo p.help
            quit QuitSuccess
    except UsageError:
        stderr.writeLine getCurrentExceptionMsg() 
        echo "Use --help for usage information"
        quit QuitSuccess
