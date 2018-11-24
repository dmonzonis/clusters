import argparse
import os

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def compute_distance_modulus(parallax):
    return 5 * np.log10(1. / (parallax / 1000.)) - 5


def load_isochrones(extinction='0.0'):
    """Return (logt, b, bprp) for the isochrone with the given Av"""
    filename = f'Av_{extinction}.dat'
    logt, g, bp, rp = np.loadtxt(filename, usecols=(1, 8, 9, 10), unpack=True)
    return logt, g, bp - rp


def load_star_data(filename):
    """Return the mean parallax, g and bp-rp of the stars in the data file."""
    ext = os.path.splitext(filename)[-1]
    if ext == '.fits':
        # Fits file
        with fits.open(filename) as hdul:
            data = hdul[1].data
        parallax, g, bp_rp = data['parallax'], data['g'], data['bp_rp']
    else:
        # ASCII file
        parallax, g, bp_rp = np.genfromtxt(filename, usecols=(5, 23, 24), unpack=True)
    return np.mean(parallax), g, bp_rp


def plot_isochrone(extinction, age, distance_modulus):
    """Plots an isochrone on top of the star data.

    Args:
        extinction: Degree of extinction, in string form. From '0.0' to '1.0'
        distance_modulus: Distance modulus to offset the stars vertically
        age: log of the age of the isochrone to plot
    """
    logti, gi, bprpi = load_isochrones(extinction)
    age_in_gyr_label = "{0:.5f}".format(10**age / 10**9)
    distance_modulus_label = "{0:.5f}".format(distance_modulus)
    plt.plot(bprpi[logti == age], gi[logti == age] + distance_modulus,
             'r-', label=f'{age_in_gyr_label} Gyr, Av={extinction}, dM={distance_modulus_label}')


def plot_star_data(g, bprp, name):
    """Create a scatter plot with the star data and return it.

    Args:
        data_filename: Name of the file with the star data
    """
    plt.title(name)
    plt.xlim(0, 3)
    plt.ylim(18, 4)
    plt.xlabel('BP-RP')
    plt.ylabel('G')
    plt.minorticks_on()
    plt.tight_layout()
    # Plot the data
    plt.scatter(bprp, g)


def main():
    # Define font for plots
    font = {'family': 'serif',
            'size': 10}
    rc('font', **font)

    # Create argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=argparse.FileType('r'),
                        help="star data file", required=True)
    parser.add_argument('-e', '--extinction', type=str,
                        help="degree of extinction", required=True)
    parser.add_argument('-d', '--distance', type=float,
                        help="distance modulus")
    parser.add_argument('-a', '--age', type=float,
                        help="log of the age of the isochrone to plot", required=True)
    parser.add_argument('-s', '--save', action='store_true',
                        help="save plot to file")
    parser.add_argument('-n', '--name', type=str,
                        help="name of the cluster shown on the plot")
    args = parser.parse_args()

    # Load star data
    filename = args.f.name
    name = args.name if args.name else os.path.splitext(filename)[0]

    mean_parallax, g, bprp = load_star_data(filename)
    distance_modulus = compute_distance_modulus(mean_parallax)

    # Create plot figure
    plt.figure(figsize=(12, 8))

    # Plot data and isochrone
    if not args.distance:
        plot_star_data(g, bprp, name)
        plot_isochrone(args.extinction, args.age, distance_modulus)
    else:
        plot_star_data(g, bprp, name)
        plot_isochrone(args.extinction, args.age, args.distance)

    # Add legends
    plt.legend()

    if args.save:
        plt.savefig(name + '.png')
    else:
        plt.show()


if __name__ == "__main__":
    main()
