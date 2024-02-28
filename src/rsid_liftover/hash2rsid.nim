import times
import strformat
import strutils
import sequtils
import streams
import tables
import ./parsers/hash2rsid
import ./utils
import ./snp_data
import std/hashes
import re

iterator readInputValues*(filename: string, sep: string, varid_idx: int, var_regexp: Regex, chrom_idx: int, pos_idx: int,  ref_idx: int, alt_idx: int, header: bool = false): SnpLine =
    var
        fileStream: Stream
        snp_data: SnpLine
        is_header_line = header
        varid_colidx = varid_idx

    fileStream = readInputFile(filename)
    defer: close(fileStream)

    log("INFO", fmt"Reading input file {filename}")
    if fileStream == nil:
        raise newException(IOError, "Cannot open file " & filename)
    
    var line: string
    while not fileStream.atEnd:
        line = fileStream.readLine()
        if is_header_line:
            is_header_line = false
            snp_data.line = line
            snp_data.id = -9
            snp_data.chrom = "HEADER"
        else: 
            let fields = line.split(sep)
            snp_data.line = line
            if fields.len == 1: varid_colidx = 0
            if varid_colidx != -1:
                var matches: array[4, string]
                if not fields[varid_colidx].match(var_regexp, matches):
                    log("ERROR", fmt"Error parsing variant id {fields[varid_colidx]}")
                    quit "", QuitFailure
                snp_data.chrom = matches[0]
                snp_data.pos = parseInt(matches[1])
                snp_data.ref_a = matches[2]
                snp_data.alt_a = matches[3]
            else:
                snp_data.id = -1
                snp_data.chrom = fields[chrom_idx].replace("chr", "")
                snp_data.pos = parseInt(fields[pos_idx])
                snp_data.alt_a = fields[alt_idx]
                snp_data.ref_a = fields[ref_idx]
            
        yield snp_data

proc getRsIdFromBin(snp: SnpLine, map_dir:string, target_build:string, selected_chroms: seq[string], chrom_data: var Table[string, ChromData]): (int, int) =
    var chr_to_search = @[snp.chrom]
    if snp.chrom == "-1": chr_to_search = selected_chroms
    
    let key_hash = hash(fmt"{snp.pos}{snp.ref_a}_{snp.alt_a}")
    var snp_id: int
    for c in chr_to_search:
        if not chrom_data.hasKey(c):
            var c_data: ChromData
            c_data.getChromData(map_dir, c, target_build, mode = "hash2rsid")
            chrom_data[c] = c_data
        snp_id = chrom_data[c].tab.getOrDefault(key_hash, -1)
        if snp_id != -1:
            return (chrom_data[c].dbsnp, snp_id)
    return (-1, -1)

proc main* (argv: seq[string]) =    
    var opts = parseCmdLine(argv)

    let
        out_prefix = opts.out
        sep = opts.sep
        varid_colidx = parseInt(opts.varid_column)
        chrom_colidx = parseInt(opts.chrom_column)
        pos_colidx = parseInt(opts.pos_column)
        alt_colidx = parseInt(opts.alt_column)
        ref_colidx = parseInt(opts.ref_column)
        target_build = opts.build
        map_dir = opts.map_dir
        header = opts.header
        skip_missing = opts.no_missing
        tool_prefix = fmt"hash2rsid_{target_build}"
        input_files = opts.intables

    if varid_colidx == -1 and any(@[chrom_colidx, pos_colidx, alt_colidx, ref_colidx], proc(x: int): bool = x == -1):
        log("ERROR", "Error: all chrom, pos, ref and alt columns must be specified when variant id column is not present")
        quit "", QuitFailure
    if varid_colidx != -1 and any(@[chrom_colidx, pos_colidx, alt_colidx, ref_colidx], proc(x: int): bool = x != -1):
        log("WARN", fmt"Variant id column is set to {varid_colidx}, chrom, pos, ref and alt columns will be ignored")

    # Compile regexp from opts
    let vid_regexp = re(opts.parse_exp)

    # Performs initial checks and return list of selected chromosomes
    let selected_chroms = initialChecks(opts.chrom, input_files) 

    #Set header for output
    var header_line = "hash2rsid_dbSNPv\thash2rsid_rsid"

    var 
        n: int
        w: int
        file_t0: float
        t0 = cpuTime()
        interval = 5000

    var hash2rsid_data: Table[string, ChromData]

    for df in input_files:
        let out_stream = getOutputStream(out_prefix, df, tool_prefix)
        file_t0 = cpuTime()
        n = 0
        w = 0
        for s in readInputValues(df, sep, varid_colidx, vid_regexp, chrom_colidx, pos_colidx, ref_colidx, alt_colidx, header=header):
            n += 1
            let (log_step, msg) = progress_counter(n, interval, file_t0)    
            if log_step: log("INFO", msg)
            if s.chrom == "HEADER":
                out_stream.writeLine(&"{s.line}\t{header_line}")
                w += 1
                continue
            
            if s.chrom != "-1" and s.chrom notin selected_chroms: continue #Only process for selected cromosomes
            
            let (dbsnp_v, target_rsid) = getRsIdFromBin(s, map_dir, target_build, selected_chroms, hash2rsid_data)
            if target_rsid == -1 and skip_missing: continue
            out_stream.writeLine(&"{s.line}\t{dbsnp_v}\trs{target_rsid}")
            w += 1
        close(out_stream)
        log("INFO", fmt"Processed {n} lines in {df}, wrote {w} lines to output")

    log("INFO", fmt"Total time: {elapsed_time(t0)}")
