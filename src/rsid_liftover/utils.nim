import times
import strutils
import strformat
import system
import math
import streams
import zip/gzipfiles
import os
import sequtils
import tables

const
    STDCHROMS* = concat(map(to_seq(1..22), proc(x: int): string = $x), @["X", "Y"])
    MAP_RSID_COLIDX* = 1
    MAP_POS_COLIDX* = 2
    MAP_REF_COLIDX* = 3
    MAP_ALT_COLIDX* = 4
    HEADER_ID* = -99

type ChromData* = object
    liftover_tab*: Table[int, int]
    chrom*: string
    #ref_a*: seq[string]
    #alt_a*: seq[string]

type AlleleData* = object
    alleles_tab*: Table[string, int] #position, (rsid, ref allele, alt allele)
    chrom*: string

#TODO Extend the above structure to include ALT alleles to manage VARID to RSID conversion

type SnpLine* = tuple
    id: int
    chrom: string
    pos: int
    ref_a: string
    alt_a: string
    line: string

proc log* (level: string = "INFO", msg: string) {.discardable.} = 
    let t = now()
    let f = initTimeFormat("HH:mm:ss-fff")
    var time_string = format(t, f)
    let log_msg = fmt"[{time_string}] - {level}: {msg}"
    stderr.write(log_msg & "\n")

proc elapsed_time* (start_time: float): string =
    let interval = cpuTime() - start_time
    let s = floor(interval)
    let m = floor(((interval - s) * 1000))
    let time_interval = initDuration(seconds = int(s), milliseconds = int(m))
    result = $time_interval

proc progress_counter* (n:int, interval: var int, t0: var float): (bool, string) {.discardable.} =
    result = (false, "")

    if floorMod(n, interval) == 0:
        result = (true, fmt"{n} variants processed. Last batch: {interval} in {elapsed_time(t0)}")
        t0 = cpuTime()

    case n
        of 50000: interval = 25000
        of 150000: interval = 50000
        of 500000: interval = 100000
        of 1000000: interval = 500000
        else: discard

proc readInputFile*(filename: string): Stream =
    result = newFileStream(filename, fmRead)

    if filename.endsWith(".gz"):
        result = newGzFileStream(filename, fmRead)

proc saveChromDataToFile*(chromData: ChromData, filename: string) =
    let fileStream = newFileStream(filename, fmWrite)
    if fileStream.isNil:
        raise newException(IOError, "Unable to open file " & filename)
    defer: fileStream.close()
    fileStream.write(chromData.chrom.len.uint8)
    fileStream.write(chromData.chrom)
    fileStream.write(chromData.liftover_tab.len.uint32)
    for x in chromData.liftover_tab.keys:
        fileStream.write(x.uint32)
    for x in chromData.liftover_tab.values:
        fileStream.write(x.uint32)

proc saveAlleleDataToFile*(alleleData: AlleleData, filename: string) =
    var n_skipped = 0
    let fileStream = newFileStream(filename, fmWrite)
    if fileStream.isNil:
        raise newException(IOError, "Unable to open file " & filename)
    defer: fileStream.close()
    fileStream.write(alleleData.chrom.len.uint8)
    fileStream.write(alleleData.chrom)
    fileStream.write(alleleData.alleles_tab.len.uint32)
    for k, v in alleleData.alleles_tab.pairs:
        if k.len > 255:
            n_skipped += 1
            log("WARN", fmt"Key {k} too long, this will not be stored")
            continue
        fileStream.write(k.len.uint8)
        fileStream.write(k)
        fileStream.write(v.uint32)
    log("INFO", fmt"Skipped {n_skipped} keys")       

proc loadChromDataFromFile*(filename: string): ChromData =
    var
        tab_len = 0'u32
        placeholder = 0'u32
    
    var f = newFileStream(filename, fmRead)
    if f == nil:
        raise newException(IOError, "could not open file:" & filename)

    # Read chrom string
    var chrom: uint8 = 0
    discard f.readData(chrom.addr, sizeof(chrom))
    result.chrom = newString(chrom)
    discard f.readData(result.chrom[0].addr, chrom.int)

    # Read rsids and pos sequences
    discard f.readData(tab_len.addr, tab_len.sizeof)
    var rsids = newSeq[int](tab_len)
    if tab_len > 0'u32:
        for x in 0..tab_len.int-1:
            discard f.readData(
                rsids[x].addr, sizeof(placeholder))
            result.liftover_tab[rsids[x]] = 0
        for x in 0..tab_len.int-1:
            discard f.readData(
                result.liftover_tab[rsids[x]].addr, sizeof(placeholder))

proc loadAlleleDataFromFile*(filename: string): AlleleData =
    var tab_len = 0'u32

    var f = newFileStream(filename, fmRead)
    if f == nil:
        raise newException(IOError, "could not open file:" & filename)

    # Read chrom string
    var chrom: uint8 = 0
    discard f.readData(chrom.addr, sizeof(chrom))
    result.chrom = newString(chrom)
    discard f.readData(result.chrom[0].addr, chrom.int)

    # Read rsids and pos sequences
    discard f.readData(tab_len.addr, tab_len.sizeof)
    echo fmt"Reading {tab_len} positions"
    if tab_len > 0'u32:
        for x in 0..tab_len.int-1:
            if floorMod(x, 500000) == 0: echo "500k processed"
            var key_len: uint8 = 0
            discard f.readData(
                key_len.addr, sizeof(key_len))
            var key_string = newString(key_len)
            discard f.readData(key_string[0].addr, key_len.int)
            var rsid: uint32 = 0
            discard f.readData(
                rsid.addr, sizeof(rsid))
            result.alleles_tab[key_string] = rsid.int
            #echo fmt"Loaded snp at position {positions[x]}: {rsid.int}, {ref_string}, {alt_string}"

proc readAlleleMap(rsid_file: string): AlleleData =
    var
        line: string
    let fileStream = readInputFile(rsid_file)
    defer: fileStream.close()

    while not fileStream.atEnd:
        line = fileStream.readLine()
        if line.startsWith("#"): continue
        let
            fields = line.split("\t")
            position = fields[MAP_POS_COLIDX]
            ref_allele = fields[MAP_REF_COLIDX] 
            alt_alleles = fields[MAP_ALT_COLIDX].split(",")
        for a in alt_alleles:
            result.alleles_tab[fmt"{position}_{ref_allele}_{a}"] = parseInt(fields[MAP_RSID_COLIDX])

proc readMap(rsid_file: string): ChromData =
    var
        line: string
    let fileStream = readInputFile(rsid_file)
    defer: fileStream.close()

    while not fileStream.atEnd:
        line = fileStream.readLine()
        if line.startsWith("#"): continue
        let
            fields = line.split("\t")
        result.liftover_tab[parseInt(fields[MAP_RSID_COLIDX])] = parseInt(fields[MAP_POS_COLIDX])

proc singleChromData*(fileprefix: string, expected_chrom: string, ignore_binaries: bool = false): ChromData =
    var t0 = cpuTime()
    let
        rsid_file_tsv = fmt"{fileprefix}.tsv.gz"
        rsid_file_bin = fmt"{fileprefix}.liftover.bin"
    if fileExists(rsid_file_bin) and not ignore_binaries:
        log("INFO", fmt"Loading from binary file {rsid_file_bin}")
        result = loadChromDataFromFile(rsid_file_bin)
        if result.chrom != expected_chrom:
            log("ERROR", fmt"{result.chrom} found in binary, {expected_chrom} was expected. Data may be corrupted")
            quit "", QuitFailure
    else: 
        log("INFO", fmt"Loading from TSV file {rsid_file_tsv}")
        if not fileExists(rsid_file_tsv):
            log("ERROR", fmt"Unable to find .tsv.gz or .bin for {fileprefix}")
            quit "", QuitFailure
        result = readMap(rsid_file_tsv)
    result.chrom = expected_chrom

    log("INFO", fmt"Data with {result.liftover_tab.len} elements loaded in {elapsed_time(t0)}")

proc singleAlleleData*(fileprefix: string, expected_chrom: string, ignore_binaries: bool = false): AlleleData =
    var t0 = cpuTime()
    let
        rsid_file_tsv = fmt"{fileprefix}.tsv.gz"
        rsid_file_bin = fmt"{fileprefix}.alleles.bin"
    if fileExists(rsid_file_bin) and not ignore_binaries:
        log("INFO", fmt"Loading from binary file {rsid_file_bin}")
        result = loadAlleleDataFromFile(rsid_file_bin)
        if result.chrom != expected_chrom:
            log("ERROR", fmt"{result.chrom} found in binary, {expected_chrom} was expected. Data may be corrupted")
            quit "", QuitFailure
    else: 
        log("INFO", fmt"Loading from TSV file {rsid_file_tsv}")
        if not fileExists(rsid_file_tsv):
            log("ERROR", fmt"Unable to find .tsv.gz or .bin for {fileprefix}")
            quit "", QuitFailure
        result = readAlleleMap(rsid_file_tsv)
    result.chrom = expected_chrom

    log("INFO", fmt"Data with {result.alleles_tab.len} elements loaded in {elapsed_time(t0)}")

proc loadChromData*(map_dir: string, target_build: string, dbsnp_version: string, selected_chroms: seq[string], ignore_binaries: bool = false): Table[string, ChromData] =
    var chrom_list = selected_chroms
    if selected_chroms[0] == "-1": chrom_list = STDCHROMS
    for chrom in chrom_list:
        var t0 = cpuTime()
        let
            file_prefix = fmt"{map_dir}/{target_build}_dbSNP{dbsnp_version}.chr{chrom}"
        result[chrom] = singleChromData(file_prefix, chrom)
        log("INFO", fmt"Chromosome {chrom} loaded in {elapsed_time(t0)}")

proc binarySearch*(arr: seq[int], target: int): int =
    var left = 0
    var right = arr.len - 1

    while left <= right:
        let mid = left + (right - left) div 2
        if arr[mid] == target:
            return mid
        elif arr[mid] < target:
            left = mid + 1
        else:
            right = mid - 1

    return -1

iterator readInputValues*(filename: string, separator: string, rsid_idx: int, chrom_idx: int, pos_idx: int = MAP_POS_COLIDX, alt_idx: int = MAP_ALT_COLIDX, ref_idx: int = MAP_REF_COLIDX, header: bool = false): SnpLine =
    var
        fileStream: Stream
        snp_data: SnpLine
        is_header_line = header

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
            snp_data.id = -9
            snp_data.chrom = "-1"
        else: 
            let fields = line.split(separator)
            snp_data.line = line
            snp_data.id = parseInt(fields[rsid_idx].replace("rs", ""))
            snp_data.chrom = "-1"
            if chrom_idx != -1 and fields.len > 1:
                snp_data.chrom = fields[chrom_idx].replace("chr", "")
            if pos_idx != -1 and fields.len > 1:
                snp_data.pos = parseInt(fields[pos_idx])
            if alt_idx != -1 and fields.len > 1:
                snp_data.alt_a = fields[alt_idx]
            if ref_idx != -1 and fields.len > 1:
                snp_data.ref_a = fields[ref_idx]
            
        yield snp_data