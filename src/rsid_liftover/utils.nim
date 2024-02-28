import times
import strutils
import sequtils
import strformat
import system
import math
import streams
import zip/gzipfiles
import os
from ./constants import STDCHROMS

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
        of 25000: interval = 10000
        of 50000: interval = 25000
        of 150000: interval = 50000
        of 500000: interval = 100000
        of 1000000: interval = 500000
        else: discard

proc getOutputStream*(outpath: string, infile_path: string, tool_prefix: string): Stream =
    if outpath.len > 0:
        var fileprefix = extractFilename(infile_path)
        fileprefix = changeFileExt(fileprefix, "")
        var outfile: string
        if dirExists(outpath):
            outfile = fmt"{outpath}/{tool_prefix}-{fileprefix}.tsv"
        else:
            outfile = fmt"{outpath}-{tool_prefix}-{fileprefix}.tsv"
        result = newFileStream(outfile, fmWrite)
        if isNil(result):
            raise newException(IOError, "Cannot open file " & outfile)
    else:
        result = newFileStream(stdout)

proc readInputFile*(filename: string): Stream =
    result = newFileStream(filename, fmRead)

    if filename.endsWith(".gz"):
        result = newGzFileStream(filename, fmRead)

proc initialChecks*(chromosomes: string, input_files: seq[string]): seq[string] =
    #Check if the input file exists
    if input_files.len == 0:
        log("ERROR", "No input files provided")
        quit "", QuitFailure
    for f in input_files:
        if not fileExists(f):
            log("ERROR", fmt"Input file '{f}' does not exist")
            quit "", QuitFailure
    
    result = chromosomes.split(",").map(proc(x: string): string = x.replace("chr", ""))
    if result.len > 1 and any(result, proc(x: string): bool = x == "-1"):
        log("ERROR", "Cannot mix -1 with other chromosomes")
        quit "", QuitFailure
    if result[0] == "-1": result = STDCHROMS
    log("INFO", fmt"Selected chromosomes: {result}")