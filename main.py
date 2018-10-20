import os
import pickle
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astroquery.gaia import Gaia
import numpy as np
import matplotlib.pyplot as plt


DATA_FOLDER = os.path.join(os.path.dirname(__file__), 'gaia_data/')


def generate_query(ra_centre, dec_centre, radius):
    return f"""SELECT ra, dec, source_id, l, b, parallax, parallax_error,
    pmra, pmra_error, pmdec, pmdec_error, ra_dec_corr, ra_parallax_corr,
    ra_pmra_corr, ra_pmdec_corr, dec_parallax_corr, dec_pmra_corr,
    dec_pmdec_corr, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr,
    phot_g_n_obs, phot_g_mean_mag as g, bp_rp
    FROM gaiadr2.gaia_source
    WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra_centre},{dec_centre},{radius}))=1
    AND phot_g_mean_mag<17;"""


def get_gaia_data(query):
    """Launch a request of the query from the Gaia servers and return the results."""
    job = Gaia.launch_job_async(query)
    return job.get_results()


def make_plots(cluster_data, data_folder=None, cluster_name='Cluster'):
    """Makes various plots with the given data, and saves to file if necessary."""
    # Initialize figure
    plt.figure(0, figsize=(8, 8))
    plt.clf()  # Necessary so the plots don't accumulate on top of each other

    # RA vs DEC plot on the top left corner
    plt.subplot(221)
    plt.title(cluster_name)
    plt.plot(cluster_data['ra'], cluster_data['dec'], 'k.')
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.xlim(min(cluster_data['ra']), max(cluster_data['ra']))
    plt.ylim(min(cluster_data['dec']), max(cluster_data['dec']))

    # BP/RP vs magnitude on the top right corner
    plt.subplot(222)
    plt.plot(cluster_data['bp_rp'], cluster_data['g'], 'k.')
    plt.xlabel('BP-RP')
    plt.ylabel('G')
    plt.xlim(-1, 3)
    plt.ylim(18, 6)

    # pmra vs pmdec on the bottom left corner
    plt.subplot(223)
    plt.plot(cluster_data['pmra'], cluster_data['pmdec'], 'k.')
    plt.xlabel('pmra')
    plt.ylabel('pmdec')

    # magnitude vs parallax on the bottom right corner
    plt.subplot(224)
    plt.plot(cluster_data['g'], cluster_data['parallax'], 'k.')
    plt.xlabel('G')
    plt.ylabel('parallax')

    plt.tight_layout()

    if data_folder is not None:
        plt.savefig(data_folder + cluster_name + '.png')
    else:
        plt.show()


def get_unverified_cluster_data(identified_filename,
                                verified_filename,
                                with_plots=True):
    """Get the unverified clusters' data, and save it to file in FITS format.

    Unverified clusters are the ones in the identified dataset that aren't in the verified
    dataset.
    The filename of the output file will be the cluster as found in the identified clusters
    dataset.
    If the with_plots option is set to True, an image file with some relevant plots will also be
    saved in the same directory.
    """
    # Data from the Sampedro file
    id_hdul = fits.open(identified_filename)  # File with identified cluster data
    id_ra = id_hdul[1].data['RA_ICRS']
    id_dec = id_hdul[1].data['DE_ICRS']
    id_name = id_hdul[1].data['Cluster']

    # Data from the verified cluster file
    ver_hdul = fits.open(verified_filename)  # File with verified cluster data
    ver_name = ver_hdul[1].data['cluster']

    ver_set = set(ver_name)  # Set of verified cluster names
    unver_set = {c for c in id_name
                 if c not in ver_set and c != "Melotte_111"}  # Melotte_111 is a well known cluster

    # Save unverified cluster names to a text file for later use
    with open('unverified_cluster_list.txt', 'w') as unver_file:
        for name in unver_set:
            unver_file.write(name + '\n')

    # Get data and plots for all the unverified clusters
    for name in unver_set:
        print(f"Getting data for cluster {name}")

        # Find centre and size of the field corresponding to the cluster
        ra_centre = np.mean([max(id_ra[id_name == name]), min(id_ra[id_name == name])])
        dec_centre = np.mean([max(id_dec[id_name == name]), min(id_dec[id_name == name])])
        radius = (max(id_dec[id_name == name]) - min(id_dec[id_name == name])) / 2.

        query = generate_query(ra_centre, dec_centre, radius)
        results = get_gaia_data(query)

        data_folder = os.path.join(dirname, 'gaia_data/')
        if not os.path.exists(data_folder):
            os.makedirs(data_folder)
        results.write(data_folder + name + '.fits', overwrite=True)
        make_plots(results, data_folder=data_folder, cluster_name=name)

    id_hdul.close()
    ver_hdul.close()


def main():
    get_unverified_cluster_data('sampedro_stars.fits', 'verified.fits')


if __name__ == "__main__":
    main()
