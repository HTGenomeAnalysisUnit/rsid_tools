import times
import strutils
import strformat
import system
import math
import streams
import zip/gzipfiles

type ChromData* = object
    rsids*: seq[int]
    chrom*: string
    pos*: seq[int]
    #ref_a*: seq[string]
    #alt_a*: seq[string]

type SnpLine* = tuple
    id: int
    chrom: string
    line: string

proc saveChromDataToFile*(chromData: ChromData, filename: string) =
    let fileStream = newFileStream(filename, fmWrite)
    if fileStream.isNil:
        raise newException(IOError, "Unable to open file " & filename)
    defer: fileStream.close()
    fileStream.write(chromData.chrom.len.uint8)
    fileStream.write(chromData.chrom)
    fileStream.write(chromData.rsids.len.uint32)
    fileStream.write(chromData.pos.len.uint32)
    for x in chromData.rsids:
        fileStream.write(x.uint32)
    for x in chromData.pos:
        fileStream.write(x.uint32)

proc loadChromDataFromFile*(filename: string): ChromData =
    var
        rsids = 0'u32
        pos = 0'u32
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
    discard f.readData(rsids.addr, rsids.sizeof)
    discard f.readData(pos.addr, pos.sizeof)
    result.rsids = newSeq[int](rsids)
    result.pos = newSeq[int](pos)
    if rsids > 0'u32:
        for x in 0..rsids.int-1:
            discard f.readData(
                result.rsids[x].addr, sizeof(placeholder))
    echo "finished rsids"
    if pos > 0'u32:
        for x in 0..pos.int-1:
            discard f.readData(
                result.pos[x].addr, sizeof(placeholder))

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

proc readInputFile*(filename: string): Stream =
    result = newFileStream(filename, fmRead)

    if filename.endsWith(".gz"):
        result = newGzFileStream(filename, fmRead)

iterator readInputValues*(filename: string, separator: string, rsid_idx: int, chrom_idx: int, header: bool): SnpLine =
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
        yield snp_data