import argparse
import os

from astropy.io import fits
import numpy as np


def load_star_data(filename):
    ext = os.path.splitext(filename)[-1]
    is_fits = ext == '.fits'
    if is_fits:
        # Fits file
        with fits.open(filename) as hdul:
            data = hdul[1].data
    else:
        # ASCII file
        data = np.genfromtxt(filename, unpack=True)
    return data, is_fits


def summary(data, is_fits):
    if is_fits:
        ra_centre = np.mean(data['ra'])
        dec_centre = np.mean(data['dec'])
        parallax = np.mean(data['parallax'])
        members = len(data['ra'])
    else:
        ra_centre = np.mean(data[0])
        dec_centre = np.mean(data[1])
        parallax = np.mean(data[5])
        members = len(data[0])

    print(f"Center: ({ra_centre}, {dec_centre})")
    print(f"Parallax: {parallax}")
    print(f"Members: {members}")


def main():
    # Create argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=argparse.FileType('r'),
                        help="star data file", required=True)
    args = parser.parse_args()

    # Load star data
    filename = args.f.name
    data, is_fits = load_star_data(filename)
    summary(data, is_fits)


if __name__ == "__main__":
    main()
