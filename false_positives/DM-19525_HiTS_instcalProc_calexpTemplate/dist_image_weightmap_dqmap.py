# Prints commands to link the fits files on the argument list into
# instcals, dqmasks, wtmaps subdirs.

import astropy.io.fits as pf
import sys
import os

def main():
    for F in sys.argv[1:]:
        with pf.open(F) as hdul:
            p = os.path.abspath(F)
            d, fname = os.path.split(p)
            lastdir = os.path.basename(d)

            x = hdul[0].header['PRODTYPE']

            if x == "image":
                target = "instcals"
            elif x == "dqmask":
                target = "dqmasks"
            elif x == "wtmap":
                target = "wtmaps"
            else:
                raise ValueError("Unknown prodtype identifier")

            print (f"ln -s {p} ../{target}/")

if __name__ == "__main__":
    main()
