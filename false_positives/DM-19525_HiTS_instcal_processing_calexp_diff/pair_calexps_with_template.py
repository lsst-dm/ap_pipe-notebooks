import astropy.io.fits as pf
import sys
import os

# '15A field': ('14A field', visit chosen as template)
FIELD_PAIRS = { 'Blind15A_26': ('Blind14_04', 289444),
                'Blind15A_42': ('Blind14A_09', 289449),
                'Blind15A_40': ('Blind14A_10', 289450)}

def main():
    for F in sys.argv[1:]:
        with pf.open(F) as hdul:
            p = os.path.abspath(F)
            d, fname = os.path.split(p)
            lastdir = os.path.basename(d)

            x = hdul[0].header['OBJECT']
            if x not in FIELD_PAIRS.keys():
                raise ValueError("Unknown field")

            v15 = hdul[0].header['EXPNUM']
            print(v15,"   ",FIELD_PAIRS[x][1])


if __name__ == "__main__":
    main()
