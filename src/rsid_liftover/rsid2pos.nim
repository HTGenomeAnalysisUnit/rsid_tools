import times
import strformat
import strutils
import sequtils
import streams
import tables
import os
import ./parsers/rsid2pos
import ./utils
import ./snp_data
import std/hashes
import ./constants

iterator readInputValues(filename: string, sep: string, rsid_idx: int, chrom_idx: int, header: bool = false): SnpLine =
    var
        fileStream: Stream
        snp_data: SnpLine
        is_header_line = header
        rsid_colidx = rsid_idx
        chrom_colidx = chrom_idx

    log("INFO", fmt"Reading input file {filename}")
    fileStream = readInputFile(filename)
    defer: close(fileStream)

    if fileStream == nil:
        raise newException(IOError, "Cannot open file " & filename)
    
    var line: string
    while not fileStream.atEnd:
        line = fileStream.readLine()
        if is_header_line:
            is_header_line = false
            snp_data.line = line
            snp_data.id = 0
            snp_data.chrom = "HEADER"
        else:
            let fields = line.split(sep)
            if fields.len == 1:
                rsid_colidx = 0
                chrom_colidx = -1
            snp_data.line = line
            try:
                snp_data.id = parseInt(fields[rsid_colidx].replace("rs", ""))
                snp_data.chrom = "-1"
                if chrom_colidx != -1:
                    snp_data.chrom = fields[chrom_colidx].replace("chr", "")
            except IndexDefect:
                log("ERROR", fmt"Error parsing line {line}. One of rsid index {rsid_colidx} or chrom index {chrom_colidx} is out of bounds")
                quit "", QuitFailure
            except ValueError:
                log("ERROR", fmt"Error parsing line {line}. The rsid value {fields[rsid_colidx]} is not a proper rsid or integer")
                quit "", QuitFailure
        yield snp_data

proc liftOverFromBin(snp: SnpLine, map_dir:string, target_build:string, selected_chroms: seq[string], chrom_data: var Table[string, ChromData]): (int, string, int) =
    var chr_to_search = @[snp.chrom]
    if snp.chrom == "-1": chr_to_search = selected_chroms
    
    var snp_pos: int
    for c in chr_to_search:
        if not chrom_data.hasKey(c):
            var c_data: ChromData
            c_data.getChromData(map_dir, c, target_build, mode = "rsid2pos")
            chrom_data[c] = c_data
        snp_pos = chrom_data[c].tab.getOrDefault(snp.id, -1)
        if snp_pos != -1:
            return (chrom_data[c].dbsnp, c, snp_pos)
    return (-1, "-1", -1)

# proc getRsIdFromBin(snp: SnpLine, map_dir:string, target_build:string, dbsnp_version: string, selected_chroms: seq[string], allele_data: var Table[string, AlleleData]): int =
#     var chr_to_search = @[snp.chrom]
#     if snp.chrom == "-1" and selected_chroms[0] == "-1": chr_to_search = STDCHROMS
#     if selected_chroms[0] != "-1": chr_to_search = selected_chroms
#     let key_hash = hash(fmt"{snp.pos}{snp.ref_a}_{snp.alt_a}")
#     for c in chr_to_search:
#         let file_prefix = fmt"{map_dir}/{target_build}_dbSNP{dbsnp_version}.chr{c}"
#         if not allele_data.hasKey(c):
#             allele_data[c] = singleAlleleData(file_prefix, c)
#         result = allele_data[c].alleles_tab.getOrDefault(key_hash, -1)

proc main* (argv: seq[string]) =    
    var opts = parseCmdLine(argv)

    let
        out_prefix = opts.out
        sep = opts.sep
        rsid_colidx = parseInt(opts.rsid_column)
        chrom_colidx = parseInt(opts.chrom_column)
        target_build = opts.target
        map_dir = opts.map_dir
        header = opts.header
        skip_missing = opts.no_missing
        tool_prefix = fmt"rsid2pos_{target_build}"

    var selected_chroms = opts.chrom.split(",").map(proc(x: string): string = x.replace("chr", ""))
    if selected_chroms.len > 1 and any(selected_chroms, proc(x: string): bool = x == "-1"):
        log("ERROR", "Cannot mix -1 with other chromosomes")
        quit "", QuitFailure
    if selected_chroms[0] == "-1": selected_chroms = STDCHROMS
    log("INFO", fmt"Selected chromosomes: {selected_chroms}")

    #Set header for output
    var header_line = &"rsid2pos_dbSNPv\trsid2pos_chrom_{target_build}\trsid2pos_pos_{target_build}"

    #Check if the input file exists
    if opts.intables.len == 0:
        log("ERROR", "No input files provided")
        quit "", QuitFailure
    for f in opts.intables:
        if not fileExists(f):
            log("ERROR", fmt"Input file '{f}' does not exist")
            quit "", QuitFailure

    var 
        n: int
        w: int
        file_t0: float
        t0 = cpuTime()
        interval = 5000

    var rsid2pos_data: Table[string, ChromData]

    for df in opts.intables:
        let out_stream = getOutputStream(out_prefix, df, tool_prefix)
        file_t0 = cpuTime()
        n = 0
        w = 0
        for s in readInputValues(df, sep, rsid_colidx, chrom_colidx, header=header):
            n += 1
            let (log_step, msg) = progress_counter(n, interval, file_t0)    
            if log_step: log("INFO", msg)
            if s.chrom == "HEADER":
                out_stream.writeLine(&"{s.line}\t{header_line}")
                continue
            
            if s.chrom != "-1" and s.chrom notin selected_chroms: continue #Only process for selected cromosomes
            
            let (dbsnp_v, target_chrom, target_pos) = liftOverFromBin(s, map_dir, target_build, selected_chroms, rsid2pos_data)
            #echo fmt"{s.id} {s.chrom} resulting in dbSNP{dbsnp_v}: chrom {target_chrom}; pos {target_pos}"
            if target_pos == -1 and skip_missing: continue
            out_stream.writeLine(&"{s.line}\t{dbsnp_v}\t{target_chrom}\t{target_pos}")
            w += 1
            # else:
            #     let target_rsid = getRsIdFromBin(s, map_dir, target_build, dbsnp_v, chrom, allele_data)
            #     if target_rsid == -1 and skip_missing: continue
            #     out_stream.writeLine(&"{s.line}\t{target_rsid}")
            #(target_chrom, target_pos) = liftOver(s, chrom_data, chrom)
        close(out_stream)
        log("INFO", fmt"Processed {n} lines in {df}, wrote {w} lines to output")

    log("INFO", fmt"Total time: {elapsed_time(t0)}")

# when isMainModule:
#     main()