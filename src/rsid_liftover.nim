import times
import strformat
import strutils
import sequtils
import streams
import tables
from os import fileExists, commandLineParams
import rsid_liftover/arg_parse
import rsid_liftover/utils

const 
    VERSION = "0.1"
    STDCHROMS = concat(map(to_seq(1..22), proc(x: int): string = $x), @["X", "Y"])
    MAP_RSID_COLIDX = 1
    MAP_POS_COLIDX = 2
    HEADER_ID = -99

proc getStream(filename: string): Stream =
    if filename.len > 0:
        result = newFileStream(filename, fmWrite)
        if isNil(result):
            raise newException(IOError, "Cannot open file " & filename)
    else:
        result = newFileStream(stdout)

proc readMap(rsid_file: string): ChromData =
    var line: string
    let fileStream = readInputFile(rsid_file)
    defer: fileStream.close()

    while not fileStream.atEnd:
        line = fileStream.readLine()
        if line.startsWith("#"): continue
        let
            fields = line.split("\t")
        
        result.rsids.add(parseInt(fields[MAP_RSID_COLIDX]))
        result.pos.add(parseInt(fields[MAP_POS_COLIDX]))

proc loadChromData(map_dir: string, target_build: string, dbsnp_version: string, selected_chroms: seq[string], ignore_binaries: bool = false): Table[string, ChromData] =
    var chrom_list = selected_chroms
    if selected_chroms[0] == "-1": chrom_list = STDCHROMS
    for chrom in chrom_list:
        var t0 = cpuTime()
        let
            rsid_file_tsv = fmt"{map_dir}/{target_build}_dbSNP{dbsnp_version}.chr{chrom}.tsv.gz"
            rsid_file_bin = fmt"{map_dir}/{target_build}_dbSNP{dbsnp_version}.chr{chrom}.bin"
        if fileExists(rsid_file_bin) and not ignore_binaries:
            log("INFO", fmt"Loading chromosome {chrom} data from binary file")
            result[chrom] = loadChromDataFromFile(rsid_file_bin)
            if result[chrom].chrom != chrom:
                log("ERROR", fmt"Chromosome {chrom} data is corrupted")
                quit "", QuitFailure
            log("INFO", fmt"Chromosome {chrom} with {result[chrom].rsids.len} elements loaded in {elapsed_time(t0)}")
            continue
        log("INFO", fmt"Loading chromosome {chrom} data from TSV file")
        if not fileExists(rsid_file_tsv):
            log("ERROR", fmt"File {rsid_file_tsv} does not exist")
            quit "", QuitFailure
        var chrom_data = readMap(rsid_file_tsv)
        
        chrom_data.chrom = chrom
        result[chrom] = chrom_data
        log("INFO", fmt"Chromosome {chrom} loaded in {elapsed_time(t0)}")

proc liftOver(snp: SnpLine, chrom_data: Table[string, ChromData], selected_chroms: seq[string]): (string, int) =
    var chr_to_search = @[snp.chrom]
    if snp.chrom == "-1" and selected_chroms[0] == "-1": chr_to_search = STDCHROMS
    if selected_chroms[0] != "-1": chr_to_search = selected_chroms
    
    for c in chr_to_search:
        let idx = binarySearch(chrom_data[c].rsids, snp.id)
        if idx < 0:
            result = ("-1", -1)
        else:
            result = (chrom_data[c].chrom, chrom_data[c].pos[idx])
            return

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
    
    var chrom = opts.chrom.split(",").map(proc(x: string): string = x.replace("chr", ""))
    if chrom.len > 1 and any(chrom, proc(x: string): bool = x == "-1"):
        log("ERROR", "Cannot mix -1 with other chromosomes")
        quit "", QuitFailure
    log("INFO", fmt"Parsed chromosomes: {chrom}")

    if make_binaries:
        log("INFO", fmt"Making binary files for chromosome data")
        var chrom_data: Table[string, ChromData]
        chrom_data = loadChromData(map_dir, target_build, dbsnp_v, chrom, ignore_binaries=true)        

        for c in chrom_data.keys():
            let bin_file = fmt"{map_dir}/{target_build}_dbSNP{dbsnp_v}.chr{c}.bin"
            log("INFO", fmt"Saving chromosome {c} data to {bin_file}")
            saveChromDataToFile(chrom_data[c], bin_file)
        quit "", QuitSuccess

    var chrom_colidx = -1
    if opts.chrom_column != "-1": chrom_colidx = parseInt(opts.chrom_column)

    #Set header for output
    let header_line = &"chrom_{target_build}\tpos_{target_build}"

    #Set output to file or stdout and write header
    var write_to_file = false
    if out_file != "": write_to_file = true

    var out_stream = getStream(out_file)

    #Check if the input file exists
    for f in opts.indf:
        if not fileExists(f):
            log("ERROR", fmt"Input file '{f}' does not exist")
            quit "", QuitFailure

    #If chrom_colidx is -1 load all chromosome data and try to map rsID
    var 
        n: int
        file_t0: float
        t0 = cpuTime()
        interval = 1000

    var chrom_data: Table[string, ChromData]
    chrom_data = loadChromData(map_dir, target_build, dbsnp_v, chrom)

    for in_x in opts.indf:
        file_t0 = cpuTime()
        n = 0
        for s in readInputValues(in_x, sep, rsid_colidx, chrom_colidx, header):
            n += 1
            progress_counter(n, interval, file_t0)
            if s.id == HEADER_ID:
                out_stream.writeLine(&"{s.line}\t{header_line}")
                continue
            let 
                (target_chrom, target_pos) = liftOver(s, chrom_data, chrom)

            out_stream.writeLine(&"{s.line}\t{target_chrom}\t{target_pos}")
           
    log("INFO", fmt"Total time: {elapsed_time(t0)}")

when isMainModule:
    main()