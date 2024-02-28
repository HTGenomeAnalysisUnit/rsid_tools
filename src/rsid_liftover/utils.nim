import times
import strutils
import strformat
import system
import math
import streams
import zip/gzipfiles
import os

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

# iterator readInputValues*(filename: string, separator: string, rsid_idx: int, chrom_idx: int, pos_idx: int = POS_COLIDX, alt_idx: int = ALT_COLIDX, ref_idx: int = REF_COLIDX, header: bool = false): SnpLine =
#     var
#         fileStream: Stream
#         snp_data: SnpLine
#         is_header_line = header

#     fileStream = readInputFile(filename)
#     defer: close(fileStream)

#     if fileStream == nil:
#         raise newException(IOError, "Cannot open file " & filename)
    
#     var line: string
#     while not fileStream.atEnd:
#         line = fileStream.readLine()
#         if is_header_line:
#             is_header_line = false
#             snp_data.line = line
#             snp_data.id = -9
#             snp_data.chrom = "-1"
#         else: 
#             let fields = line.split(separator)
#             snp_data.line = line
#             snp_data.id = parseInt(fields[rsid_idx].replace("rs", ""))
#             snp_data.chrom = "-1"
#             if chrom_idx != -1 and fields.len > 1:
#                 snp_data.chrom = fields[chrom_idx].replace("chr", "")
#             if pos_idx != -1 and fields.len > 1:
#                 snp_data.pos = parseInt(fields[pos_idx])
#             if alt_idx != -1 and fields.len > 1:
#                 snp_data.alt_a = fields[alt_idx]
#             if ref_idx != -1 and fields.len > 1:
#                 snp_data.ref_a = fields[ref_idx]
            
#         yield snp_data