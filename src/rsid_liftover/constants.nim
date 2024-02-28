import sequtils

const
    VERSION* = "0.1.0"
    STDCHROMS* = concat(map(to_seq(1..22), proc(x: int): string = $x), @["X", "Y"])
    CHROM_COLIDX* = 0
    RSID_COLIDX* = 1
    POS_COLIDX* = 2
    REF_COLIDX* = 3
    ALT_COLIDX* = 4
    HEADER_ID* = -99