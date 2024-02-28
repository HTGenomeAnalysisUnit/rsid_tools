import argparse
import strformat
from ../utils import log
import ../constants

var p = newParser("make_bin"):
    help("make_bin: generate bin file for a dbSNP variant table")
    arg("intables", nargs = -1, help="dbSNP tables to convert to binary")
    option("-b", "--build", help="Genome build", required = true, choices = @["GRCh37", "GRCh38"])
    option("-v", "--dbsnp_version", help="dbSNP version", required = true)
    option("-o", "--out", help="Output folder", default = some("./"))
    option("-i", "--rsid_column", help="0-based index of column containing rsIDs in the input", default = some($RSID_COLIDX))
    option("-x", "--chrom_column", help="0-based index of column containing chromosome in the input", default = some($CHROM_COLIDX))
    option("-p", "--pos_column", help="0-based index of column containing position in the input", default = some($POS_COLIDX))
    option("-r", "--ref_column", help="0-based index of column containing REF allele in the input", default = some($REF_COLIDX))
    option("-a", "--alt_column", help="0-based index of column containing ALT allele(s) in the input. If multiple ALTs they must be comma-separated", default = some($ALT_COLIDX))
    run:
        log("ARG", fmt"Input file: {opts.intables.len}")
        log("ARG", fmt"Genome build: {opts.build}")
        log("ARG", fmt"dbSNP version: {opts.dbsnp_version}")
        log("ARG", fmt"Output folder: {opts.out}")
        log("ARG", fmt"rsID col index: {opts.rsid_column}")
        log("ARG", fmt"Chrom col index: {opts.chrom_column}")
        log("ARG", fmt"Pos col index: {opts.pos_column}")
        log("ARG", fmt"REF col index: {opts.ref_column}")
        log("ARG", fmt"ALT col index: {opts.alt_column}")

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
