import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def load_isochrones(extinction='0.0'):
    """Return (logt, b, bprp) for the isochrone with the given Av"""
    filename = f'Av_{extinction}.dat'
    logt, g, bp, rp = np.loadtxt(filename, usecols=(1, 8, 9, 10), unpack=True)
    return logt, g, bp - rp


def load_star_data(filename):
    """Return the data in cols 23 and 24, assumed to be g and bprp"""
    #  return np.loadtxt(filename, usecols=(23, 24), unpack=True)
    return np.genfromtxt(filename, usecols=(23, 24), unpack=True)


def plot_isochrone(extinction, distance_modulus, age):
    """Plots an isochrone on top of the star data.

    Args:
        extinction: Degree of extinction, in string form. From '0.0' to '1.0'
        distance_modulus: Distance modulus to offset the stars vertically
        age: log of the age of the isochrone to plot
    """
    logti, gi, bprpi = load_isochrones(extinction)
    # TODO: Add label
    age_in_gyr = "{0:.5f}".format(10**age / 10**9)
    plt.plot(bprpi[logti == age], gi[logti == age] + distance_modulus,
             'r-', label=f'{age_in_gyr} Gyr, Av={extinction}, dM={distance_modulus}')


def create_plot(data_filename):
    """Create a scatter plot with the star data and return it.

    Args:
        data_filename: Name of the file with the star data
    """
    # Set up plot figure
    plt.figure(figsize=(12, 8))
    # TODO: Use better plot title
    plt.title(data_filename)
    plt.xlim(0, 3)
    plt.ylim(18, 4)
    plt.xlabel('BP-RP')
    plt.ylabel('G')
    plt.minorticks_on()
    plt.tight_layout()
    # Load and plot the data
    g, bprp = load_star_data(data_filename)
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
                        help="distance modulus", required=True)
    parser.add_argument('-a', '--age', type=float,
                        help="log of the age of the isochrone to plot", required=True)
    parser.add_argument('-s', '--save', action='store_true',
                        help="save plot to file")
    args = parser.parse_args()

    # Plot data and isochrone
    create_plot(args.f.name)
    plot_isochrone(args.extinction, args.distance, args.age)

    # Add legends
    plt.legend()

    if args.save:
        plt.savefig(os.path.splitext(args.f.name)[0] + '.png')
    else:
        plt.show()


if __name__ == "__main__":
    main()
