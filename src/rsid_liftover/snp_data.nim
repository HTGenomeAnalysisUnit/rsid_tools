import times
import strutils
import strformat
import system
import streams
import zip/gzipfiles
import sequtils
import tables
import std/hashes
import ./utils
import ./constants
import sets

const MODALITIES = {
    1: "rsid2pos",
    2: "hash2rsid"
}.toTable

type ChromData* = object
    chrom*: string
    build*: string #genome build
    dbsnp*: int #dbSNP version
    mode*: int #1: rsid2pos, 2: hash2rsid
    tab*: Table[int, int] #rsid numeric id, pos

proc saveToBin*(chromData: ChromData, target_dir: string) =
    let filename = fmt"{target_dir}/{chromData.build}_{chromData.chrom}.{MODALITIES[chromData.mode]}.bin"
    log("INFO", fmt"Writing data to binary file {filename}")
    let fileStream = newFileStream(filename, fmWrite)
    if fileStream.isNil:
        raise newException(IOError, "Unable to open file " & filename)
    defer: fileStream.close()
    fileStream.write(chromData.chrom.len.uint8)
    fileStream.write(chromData.chrom)
    fileStream.write(chromData.build.len.uint8)
    fileStream.write(chromData.build)
    fileStream.write(chromData.dbsnp.uint8)
    fileStream.write(chromData.mode.uint8)
    fileStream.write(chromData.tab.len.uint32)
    for k, v in chromData.tab.pairs:
        fileStream.write(k.uint32)
        fileStream.write(v.uint32)

proc loadFromBin(c: var ChromData, filename: string) =
    log("INFO", fmt"Reading data from binary file {filename}")
    var
        tab_len = 0'u32
        u32 = 0'u32
        u8 = 0'u8
        t0 = cpuTime()
    
    var f = newFileStream(filename, fmRead)
    if f == nil:
        raise newException(IOError, "could not open file:" & filename)

    # Read chrom string
    var chrom: uint8 = 0
    discard f.readData(chrom.addr, sizeof(chrom))
    c.chrom = newString(chrom)
    discard f.readData(c.chrom[0].addr, chrom.int)

    # Read build string
    var build: uint8 = 0
    discard f.readData(build.addr, sizeof(build))
    c.build = newString(build)
    discard f.readData(c.build[0].addr, build.int)

    # Read version and mode
    discard f.readData(c.dbsnp.addr, sizeof(u8))
    discard f.readData(c.mode.addr, sizeof(u8))

    # Read keys and values for the table
    discard f.readData(tab_len.addr, tab_len.sizeof)
    var rsids = newSeq[int](tab_len)
    var positions = newSeq[int](tab_len)
    var res_tab = initTable[int, int](tab_len.int)
    #echo fmt"Reading {tab_len} elements"
    if tab_len > 0'u32:
        for x in 0..tab_len.int-1:
            discard f.readData(
                rsids[x].addr, sizeof(u32))
            discard f.readData(
                positions[x].addr, sizeof(u32))

        for pairs in zip(rsids, positions):
            let (rsid, pos) = pairs
            res_tab[rsid] = pos
        c.tab = res_tab
    log("INFO", fmt"Loaded {tab_len} elements in {elapsed_time(t0)}")

proc loadFromTsv*(filename: string, build: string, dbsnp_v, rsid_colidx, chrom_colidx, pos_colidx, ref_colidx, alt_colidx: int): (string, ChromData, ChromData) =
    log("INFO", fmt"Reading SNP file {filename}")
    var 
        line: string
        rsid_to_pos: ChromData
        hash_to_rsid: ChromData
        chrom_list: HashSet[string]
        n_skipped = 0
        t0 = cpuTime()

    rsid_to_pos.build = build
    rsid_to_pos.mode = 1
    rsid_to_pos.dbsnp = dbsnp_v
    hash_to_rsid.build = build
    hash_to_rsid.mode = 2
    hash_to_rsid.dbsnp = dbsnp_v
    let fileStream = readInputFile(filename)
    defer: fileStream.close()

    while not fileStream.atEnd:
        line = fileStream.readLine()
        if line.startsWith("#"): continue
        let
            fields = line.split("\t")
            chrom = fields[chrom_colidx].replace("chr", "")
            rsid = parseInt(fields[rsid_colidx].replace("rs", ""))
            position = parseInt(fields[pos_colidx])
            ref_allele = fields[ref_colidx] 
            alt_alleles = fields[alt_colidx].split(",")
            
        chrom_list.incl(chrom)
        if chrom_list.len > 1:
            log("ERROR", fmt"Multiple chromosomes {chrom_list} found in the file, this is not supported")
            quit "", QuitFailure

        rsid_to_pos.tab[rsid] = position
        for a in alt_alleles:
            let
                k = fmt"{position}{ref_allele}_{a}"
                k_enc = hash(k)
            #echo fmt"Encoded {k} as {k_enc}"
            if k.len > 255:
                n_skipped += 1
                log("WARN", fmt"Key {k} too long, this will not be loaded")
                continue
            if k_enc > 4294967295:
                n_skipped += 1
                log("WARN", fmt"Key {k} can not be encoded in 32bit hash, this will not be loaded")
                continue
            hash_to_rsid.tab[k_enc] = rsid
    let chrom_value = chrom_list.pop()
    rsid_to_pos.chrom = chrom_value
    hash_to_rsid.chrom = chrom_value
    log("INFO", fmt"Loaded {rsid_to_pos.tab.len} elements in {elapsed_time(t0)}")
    log("INFO", fmt"Skipped {n_skipped} SNPs in the allele table")
    return (chrom_value, rsid_to_pos, hash_to_rsid)

proc getChromData*(c: var ChromData, data_dir: string, chrom: string, build: string, mode: string) =
    let bin_file = fmt"{data_dir}/{build}_{chrom}.{mode}.bin"
    c.loadFromBin(bin_file)
    if c.chrom != chrom:
        log("ERROR", fmt"{c.chrom} found in {bin_file}, {chrom} was expected. Data may be corrupted")
        quit "", QuitFailure
    if c.build != build:
        log("ERROR", fmt"{c.build} found in {bin_file}, {build} was expected. Data may be corrupted")
        quit "", QuitFailure

proc getMultipleChrom*(data_dir: string, build: string, selected_chroms: seq[string], mode: string): Table[string, ChromData] =
    var chrom_list = selected_chroms
    if selected_chroms[0] == "-1": chrom_list = STDCHROMS
    for chrom in chrom_list:
        var c_data: ChromData
        c_data.getChromData(data_dir, chrom, build, mode)
        result[chrom] = c_data
