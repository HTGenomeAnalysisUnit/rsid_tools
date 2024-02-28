import argparse
import os
import tables
import strformat
from ./rsid_liftover/utils import log
from ./rsid_liftover/constants import VERSION
import ./rsid_liftover/rsid2pos
# import ./rsid_liftover/var2rsid
import ./rsid_liftover/make_bin

proc main*() =
  type pair = object
    f: proc(argv: seq[string])
    description: string

  var dispatcher = {
    "liftover": pair(f:rsid2pos.main, description:"get variant coordinates based on rsIDs"),
    # "annotate": pair(f:var2rsid.main, description:"annotate rsIDs for a list of variants"),
    "make_bin": pair(f:make_bin.main, description:"generate binaries files from a dbSNP table")
    }.toOrderedTable

  var args = commandLineParams()
  
  log("INFO", fmt"rsID tools v{VERSION}")

  if len(args) == 0 or args[0] in @["-h", "--help"]:
    stderr.write_line "\nCommands: "
    for k, v in dispatcher:
      echo &"  {k:<13}:   {v.description}"
    quit ""
  if len(args) > 0 and (args[0] notin dispatcher):
    echo &"unknown program '{args[0]}'"
    quit ""

  var tool_args = args[1..args.high]
  if tool_args.len == 0:
    tool_args = @["--help"]
  dispatcher[args[0]].f(tool_args)

when isMainModule:
  main()