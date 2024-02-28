import argparse
import strformat
from ../utils import log

var p = newParser("make_bin"):
    help("make_bin: generate bin file for a dbSNP variant table")
    arg("intables", nargs = -1, help="dbSNP tables to convert to binary")
    option("-s", "--sep", help="Column separator", default = some("\t"))
    flag("-e", "--header", help="Input file has header")
    option("-b", "--build", help="Genome build", required = true, choices = @["GRCh37", "GRCh38"])
    option("-o", "--out", help="Output folder or prefix. Output structure is: rsid2pos_{target_build}-{infileprefix}.tsv. If not provided output to stdout")
    option("-i", "--varid_column", help="0-based index of column containing variantID. If set all others columns will be ignored and information are parsed from id values", default = some("-1"))
    option("-e", "--parse_exp", help="Pattern for parsing chrom, pos, ref, alt from variant ID", default=some("([0-9XYM]+):([0-9]+):([ACTG]+):([ACTG]+)"))
    option("-c", "--chrom", help="Only process data for specific chromosomes. Comma-sep list accepted or -1 for all chromosomes", default = some("-1"))
    option("-x", "--chrom_column", help="0-based index of column containing chromosome in the input", default = some("-1"))
    option("-p", "--pos_column", help="0-based index of column containing position in the input", default = some("-1"))
    option("-r", "--ref_column", help="0-based index of column containing REF allele in the input", default = some("-1"))
    option("-a", "--alt_column", help="0-based index of column containing ALT allele in the input", default = some("-1"))
    option("-d", "--map_dir", help="Directory containing binaries files generated by make_bin", required = true)
    flag("-n", "--no_missing", help="Do not include in output rsIDs that can't be liftovered to the target build")
    run:
        log("ARG", fmt"Input file: {opts.intables.len}")
        log("ARG", fmt"Column separator: {opts.sep}")
        log("ARG", fmt"Header: {opts.header}")
        log("ARG", fmt"Chromosomes: {opts.chrom}")
        log("ARG", fmt"Genome build: {opts.build}")
        log("ARG", fmt"Output folder: {opts.out}")
        log("ARG", fmt"Variant ID col index: {opts.varid_column}")
        log("ARG", fmt"Parse expression: {opts.parse_exp}")
        log("ARG", fmt"Chrom col index: {opts.chrom_column}")
        log("ARG", fmt"Pos col index: {opts.pos_column}")
        log("ARG", fmt"REF col index: {opts.ref_column}")
        log("ARG", fmt"ALT col index: {opts.alt_column}")
        log("ARG", fmt"Map dir: {opts.map_dir}")

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
