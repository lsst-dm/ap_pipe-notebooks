import astropy.io.fits as pf
import sys

def main():
    for F in sys.argv[1:]:
        hdul = pf.open(F)
        print (hdul[0].header['EXPNUM']," ",F)

if __name__ == "__main__":
    main()
