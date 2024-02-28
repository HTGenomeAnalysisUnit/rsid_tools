import times
import strformat
import strutils
import os
import ./parsers/make_bin
import ./utils
import ./snp_data

proc main* (argv: seq[string]) =    
    var opts = parseCmdLine(argv)

    let
        rsid_colidx = parseInt(opts.rsid_column)
        chrom_colidx = parseInt(opts.chrom_column)
        pos_colidx = parseInt(opts.pos_column)
        ref_colidx = parseInt(opts.ref_column)
        alt_colidx = parseInt(opts.alt_column)
        dbsnp_v = parseInt(opts.dbsnp_version)
        genome_build = opts.build
        out_dir = opts.out
    
    if not dirExists(out_dir):
        echo "Output directory does not exist. Creating it..."
        createDir(out_dir)

    var t0 = cpuTime()
    log("INFO", fmt"Making binary files for SNP data from {opts.intables.len} input tables...")
    for df in opts.intables:
        let (rsid2pos, hash2rsid) = loadFromTsv(df, genome_build, dbsnp_v, rsid_colidx, chrom_colidx, pos_colidx, ref_colidx, alt_colidx)
        rsid2pos.saveToBin(out_dir)
        hash2rsid.saveToBin(out_dir)

    log("INFO", fmt"Conversion completed in {elapsed_time(t0)}")