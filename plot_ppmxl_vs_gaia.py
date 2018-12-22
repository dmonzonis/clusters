import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


DEG_TO_MAS = 3600000.  # 1 deg = 3600000 mas


def load_data():
    with fits.open('ivanov_2_ppmxl.fits') as hdul:
        ppxml_data = hdul[1].data

    gaia_data = np.genfromtxt('ivanov_2_gaiadr2.csv', skip_header=1, delimiter=',', unpack=True)
    return ppxml_data, gaia_data


def plot_positions(ppxml_data, gaia_data):
    plt.figure(0, figsize=(12, 6))
    plt.clf()
    plt.title('PPMXL vs Gaia DR2 positions')

    # PPMXL
    plt.subplot(121)
    plt.plot(ppxml_data['raj2000'], ppxml_data['dej2000'], 'k.', markersize=4)
    plt.errorbar(ppxml_data['raj2000'], ppxml_data['dej2000'],
                 xerr=ppxml_data['e_raepRA'], yerr=ppxml_data['e_deepDE'],
                 linestyle='', color='lightblue', marker='', zorder=0)
    plt.xlabel('ra (deg)')
    plt.ylabel('dec (deg)')

    # Gaia DR2 (filtering by G <= 19)
    plt.subplot(122)
    plt.plot(gaia_data[0], gaia_data[2], 'k.', markersize=4)
    # ra and dec errors are in mas, not deg
    plt.errorbar(gaia_data[0], gaia_data[2],
                 xerr=gaia_data[1] / DEG_TO_MAS, yerr=gaia_data[3] / DEG_TO_MAS,
                 linestyle='', color='lightblue', marker='', zorder=0)
    plt.xlabel('ra (deg)')
    plt.ylabel('dec (deg)')

    plt.tight_layout()

    plt.show()
    plt.close()


def plot_motions(ppxml_data, gaia_data):
    plt.figure(0, figsize=(6, 12))
    plt.clf()
    plt.title('PPMXL vs Gaia DR2 proper motions')

    # PPMXL
    plt.subplot(211)
    # In PPMXL, pmra and pmdec is given in deg/yr, so convert it to mas/yr
    plt.plot(ppxml_data['pmRA'] * DEG_TO_MAS, ppxml_data['pmDE'] * DEG_TO_MAS, 'k.', markersize=4)
    plt.errorbar(ppxml_data['pmRA'] * DEG_TO_MAS, ppxml_data['pmDE'] * DEG_TO_MAS,
                 xerr=ppxml_data['e_pmRA'] * DEG_TO_MAS, yerr=ppxml_data['e_pmDE'] * DEG_TO_MAS,
                 linestyle='', color='lightblue', marker='', zorder=0)
    plt.xlabel('pmra (mas/yr)')
    plt.ylabel('pmdec (mas/yr)')
    plt.xlim(-50, 50)
    plt.ylim(-50, 50)

    # Gaia DR2
    plt.subplot(212)
    plt.plot(gaia_data[4], gaia_data[6], 'k.', markersize=4)
    plt.errorbar(gaia_data[4], gaia_data[6],
                 xerr=gaia_data[5], yerr=gaia_data[7],
                 linestyle='', color='lightblue', marker='', zorder=0)
    plt.xlabel('pmra (mas/yr)')
    plt.ylabel('pmdec (mas/yr)')
    plt.xlim(-50, 50)
    plt.ylim(-50, 50)

    plt.tight_layout()

    plt.show()
    plt.close()


def main():
    ppxml_data, gaia_data = load_data()
    plot_positions(ppxml_data, gaia_data)
    plot_motions(ppxml_data, gaia_data)


if __name__ == "__main__":
    main()
