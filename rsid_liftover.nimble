# Package

version       = "0.1"
author        = "edoardo.giacopuzzi"
description   = "Fast liftover of SNPs based on rsID"
license       = "MIT"
srcDir        = "src"
binDir        = "bin"
bin           = @["rsid_liftover"]
skipDirs      = @["test"]


# Dependencies

requires "nim >= 1.4.8", "argparse >= 3.0.0", "zip"
