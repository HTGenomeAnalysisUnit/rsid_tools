import times
import strformat
import strutils
import sequtils
import streams
import tables
from os import fileExists, commandLineParams, extractFilename, changeFileExt, dirExists
import rsid_liftover/arg_parse
import rsid_liftover/utils

const 
    VERSION = "0.1"

proc getStream(outpath: string, in_file_path: string, targer_v: string): Stream =
    if outpath.len > 0:
        var fileprefix = extractFilename(in_file_path)
        fileprefix = changeFileExt(fileprefix, "")
        var outfile: string
        if dirExists(outpath):
            outfile = fmt"{outpath}/{fileprefix}.{targer_v}.tsv"
        else:
            outfile = fmt"{outpath}_{fileprefix}.{targer_v}.tsv"
        result = newFileStream(outfile, fmWrite)
        if isNil(result):
            raise newException(IOError, "Cannot open file " & outfile)
    else:
        result = newFileStream(stdout)

proc loadValueFromBin*(s: SnpLine, filename: string): (string, int) =
    var
        #rsids = 0'u32
        pos = 0'u32
        placeholder = 0'u32
    
    var f = newFileStream(filename, fmRead)
    if f == nil:
        raise newException(IOError, "could not open file:" & filename)

    # Read chrom string
    var chrom: uint8 = 0
    discard f.readData(chrom.addr, sizeof(chrom))
    var chrom_from_bin = newString(chrom)
    discard f.readData(chrom_from_bin[0].addr, chrom.int)
    if s.chrom != chrom_from_bin and s.chrom != "-1":
        return (s.chrom, -1)

    # Read rsids and pos sequences  
    #discard f.readData(rsids.addr, rsids.sizeof)
    discard f.readData(pos.addr, pos.sizeof)
    #var rsids_from_bin = newSeq[int](rsids)
    var pos_from_bin = newSeq[int](pos)
    # if rsids > 0'u32:
    #     for x in 0..rsids.int-1:
    #         discard f.readData(
    #             rsids_from_bin[x].addr, sizeof(placeholder))
    if pos > 0'u32:
        for x in 0..pos.int-1:
            discard f.readData(
                pos_from_bin[x].addr, sizeof(placeholder))
    echo fmt"for SNP id {s.id} found pos {pos_from_bin[s.id]} in binary file"

    result = (s.chrom, pos_from_bin[s.id])

proc liftOverFromBin(snp: SnpLine, map_dir:string, target_build:string, dbsnp_version: string, selected_chroms: seq[string], chrom_data: var Table[string, ChromData]): (string, int) =
    var chr_to_search = @[snp.chrom]
    if snp.chrom == "-1" and selected_chroms[0] == "-1": chr_to_search = STDCHROMS
    if selected_chroms[0] != "-1": chr_to_search = selected_chroms
    
    for c in chr_to_search:
        let file_prefix = fmt"{map_dir}/{target_build}_dbSNP{dbsnp_version}.chr{c}"
        if not chrom_data.hasKey(c):
            chrom_data[c] = singleChromData(file_prefix, c)

        result = (c, chrom_data[c].liftover_tab.getOrDefault(snp.id, -1))
        #result = loadValueFromBin(snp, rsid_file_bin)
        

# proc liftOver(snp: SnpLine, chrom_data: Table[string, ChromData], selected_chroms: seq[string]): (string, int) =
#     var chr_to_search = @[snp.chrom]
#     if snp.chrom == "-1" and selected_chroms[0] == "-1": chr_to_search = STDCHROMS
#     if selected_chroms[0] != "-1": chr_to_search = selected_chroms
    
#     for c in chr_to_search:
#         let idx = binarySearch(chrom_data[c].rsids, snp.id)
#         if idx < 0:
#             result = ("-1", -1)
#         else:
#             result = (chrom_data[c].chrom, chrom_data[c].pos[idx])
#             return

proc main* () =
    log("INFO", fmt"rsID liftOver v{VERSION}")

    var make_binaries = false
    var argv = commandLineParams()

    #Parse command line arguments
    if len(argv) > 0 and argv[0] == "make_chrom_binary":
        make_binaries = true

    if len(argv) == 0: 
        argv = @["--help"]
    
    var opts = parseCmdLine(argv)
    opts.logArgs()

    let
        out_file = opts.out
        sep = opts.sep
        rsid_colidx = parseInt(opts.rsid_column)
        target_build = opts.target
        map_dir = opts.map_dir
        dbsnp_v = opts.version
        header = opts.header
        skip_missing = opts.no_missing
    
    var chrom = opts.chrom.split(",").map(proc(x: string): string = x.replace("chr", ""))
    if chrom.len > 1 and any(chrom, proc(x: string): bool = x == "-1"):
        log("ERROR", "Cannot mix -1 with other chromosomes")
        quit "", QuitFailure
    log("INFO", fmt"Parsed chromosomes: {chrom}")

    if make_binaries:
        log("INFO", fmt"Making binary files for chromosome data")
        var chrom_data: ChromData     

        for c in chrom:
            log("INFO", fmt"Loading chromosome {c} data")
            let file_prefix = fmt"{map_dir}/{target_build}_dbSNP{dbsnp_v}.chr{c}"
            let bin_file = fmt"{file_prefix}.bin"
            chrom_data = singleChromData(file_prefix, c)
            log("INFO", fmt"Saving chromosome {c} data to {bin_file}")
            saveChromDataToFile(chrom_data, bin_file)
        quit "", QuitSuccess

    var chrom_colidx = -1
    if opts.chrom_column != "-1": chrom_colidx = parseInt(opts.chrom_column)

    #Set header for output
    let header_line = &"chrom_{target_build}\tpos_{target_build}"

    #Set output to file or stdout and write header
    var write_to_file = false
    if out_file != "": write_to_file = true

    #Check if the input file exists
    for f in opts.indf:
        if not fileExists(f):
            log("ERROR", fmt"Input file '{f}' does not exist")
            quit "", QuitFailure

    #If chrom_colidx is -1 load all chromosome data and try to map rsID
    var 
        n: int
        w: int
        file_t0: float
        t0 = cpuTime()
        interval = 1000

    var chrom_data: Table[string, ChromData]
    #chrom_data = loadChromData(map_dir, target_build, dbsnp_v, chrom)

    for in_x in opts.indf:
        var out_stream = getStream(out_file, in_x, dbsnp_v)
        file_t0 = cpuTime()
        n = 0
        w = 0
        for s in readInputValues(in_x, sep, rsid_colidx, chrom_colidx, header):
            n += 1
            let (log_step, msg) = progress_counter(n, interval, file_t0)    
            if log_step: log("INFO", msg)
            if s.id == HEADER_ID:
                out_stream.writeLine(&"{s.line}\t{header_line}")
                continue
            let 
                (target_chrom, target_pos) = liftOverFromBin(s, map_dir, target_build, dbsnp_v, chrom, chrom_data)
                #(target_chrom, target_pos) = liftOver(s, chrom_data, chrom)

            if target_pos == -1 and skip_missing: continue
            out_stream.writeLine(&"{s.line}\t{target_chrom}\t{target_pos}")
            w += 1
        log("INFO", fmt"Processed {n} lines in {in_x}, wrote {w} lines to output")

    log("INFO", fmt"Total time: {elapsed_time(t0)}")

when isMainModule:
    main()