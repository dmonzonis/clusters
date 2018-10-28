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

    plt.close()


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
    with fits.open(identified_filename) as hdul:
        id_data = hdul[1].data

    # Data from the verified cluster file
    with fits.open(verified_filename) as hdul:
        ver_data = hdul[1].data

    ver_set = set(ver_data['cluster'])  # Set of verified cluster names
    unver_set = {c for c in id_data['Cluster']
                 if c not in ver_set and c != "Melotte_111"}  # Melotte_111 is a well known cluster

    # Save unverified cluster names to a text file for later use
    with open('unverified_cluster_list.txt', 'w') as unver_file:
        for name in unver_set:
            unver_file.write(name + '\n')

    # Get data and plots for all the unverified clusters
    for name in unver_set:
        print(f"Getting data for cluster {name}")

        # Find centre and size of the field corresponding to the cluster
        ra_centre = np.mean([max(id_data['ra'][id_data['Cluster'] == name]),
                             min(id_data['ra'][id_data['Cluster'] == name])])
        dec_centre = np.mean([max(id_data['dec'][id_name == name]),
                              min(id_data['dec'][id_data['Cluster'] == name])])
        radius = (max(id_data['dec'][id_name == name]) -
                  min(id_data['dec'][id_data['Cluster'] == name])) / 2.

        query = generate_query(ra_centre, dec_centre, radius)
        results = get_gaia_data(query)

        if not os.path.exists(DATA_FOLDER):
            os.makedirs(DATA_FOLDER)
        results.write(DATA_FOLDER + name + '.fits', overwrite=True)
        make_plots(results, data_folder=DATA_FOLDER, cluster_name=name)


def match_cluster(cluster_name, identified_filename, min_class_matches):
    """Returns the matching of Sampedro's data with GAIA's data.

    The matches are returned in a list.
    Each match has the format (idx, d2d, d3d), where idx is the index into GAIA's data
    which corresponds to the closest object to each of the coordinates in Sampedro's data,
    d2d is the on-sky distance between them, and d3d the 3D distance.
    """
    try:
        with fits.open(DATA_FOLDER + cluster_name + '.fits') as hdul:
            cluster_data = hdul[1].data
    except FileNotFoundError:
        print("Gaia data not found. Make sure to retrieve it first.")
        return

    with fits.open(identified_filename) as hdul:
        id_data = hdul[1].data

    cluster_ra = cluster_data['ra']
    cluster_dec = cluster_data['dec']

    # Get the condition to filter the data
    # TODO: Create a good condition to filter enough stars but not too many
    condition = id_data['Cluster'] == cluster_name
    if min_class_matches == 2:
        subcondition = (id_data['ClassM1'] == 1) & (id_data['ClassM2'] == 1)
        subcondition |= (id_data['ClassM1'] == 1) & (id_data['ClassM3'] == 1)
        subcondition |= (id_data['ClassM2'] == 1) & (id_data['ClassM3'] == 1)
        condition &= subcondition
    elif min_class_matches == 3:
        condition &= id_data['ClassM1'] == 1
        condition &= id_data['ClassM2'] == 1
        condition &= id_data['ClassM3'] == 1
    else:
        raise ValueError("Can only set to match the 3 classes or at least 2")

    # Get only the data for the relevant cluster
    id_ra = id_data['RA_ICRS'][condition]
    print(len(id_ra))
    id_dec = id_data['DE_ICRS'][condition]

    gaia_star = SkyCoord(cluster_ra * u.deg, cluster_dec * u.deg, unit=(u.degree, u.degree))
    sampedro_star = SkyCoord(id_ra, id_dec, unit=(u.degree, u.degree))
    matches = sampedro_star.match_to_catalog_sky(gaia_star)

    return matches


def get_all_matches(min_class_matches=2):
    try:
        with open('unverified_cluster_list.txt') as f:
            clusters = [c.strip() for c in f.readlines()]
    except FileNotFoundError:
        print("Gaia data not found. Make sure to retrieve it first.")

    i = 1
    total = len(clusters)
    for cluster in clusters:
        print(f"Matching cluster {cluster}. Cluster {i}/{total}.")
        matches = match_cluster(cluster, "sampedro_stars.fits", min_class_matches)
        with open(DATA_FOLDER + cluster + '_match.dat', 'wb') as f:
            pickle.dump(matches, f)
        i += 1


def plot_match(cluster_name, verified_filename, data_folder=DATA_FOLDER, figsize=(12, 12)):
    with fits.open(DATA_FOLDER + cluster_name + '.fits') as hdul:
        data = hdul[1].data

    # Indices in the gaia data that correspond to
    with open(DATA_FOLDER + cluster_name + '_match.dat', 'rb') as f:
        indices = pickle.load(f)[0]

    if len(indices) == 0:
        print("No matches for cluster", cluster_name)
        with open('fail.txt', 'a') as err:
            err.write(cluster_name + '\n')
        return

    # Plot stuff
    plt.figure(figsize=figsize)
    plt.suptitle(cluster_name)

    # (ra,dec)
    plt.subplot(221)
    # All stars by Gaia
    plt.plot(data['ra'], data['dec'], '.', markersize=1, color='gray')
    # Stars also in the Sampedro catalogue
    plt.plot(data['ra'][indices], data['dec'][indices], '.', markersize=3, color='red')
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.xlim(np.nanmin(data['ra']), np.nanmax(data['ra']))
    plt.ylim(np.nanmin(data['dec']), np.nanmax(data['dec']))

    # Plot all verified cluster positions to see if one ends up in the graphic, to match
    # possible verified clusters which had different names in Sampedro's data
    with fits.open(verified_filename) as hdul:
        verified_data = hdul[1].data
    plt.plot(verified_data['ra'], verified_data['dec'], 'yX', markersize=20)
    # Add the names as labels
    for i in range(len(verified_data['cluster'])):
        xy = (verified_data['ra'][i], verified_data['dec'][i])
        plt.annotate(verified_data['cluster'][i], xy, weight='bold')

    # (pmra,pmdec)
    plt.subplot(222)
    # All stars by Gaia
    plt.plot(data['pmra'], data['pmdec'], '.', markersize=1, color='gray')
    # Stars also in the Sampedro catalogue
    plt.plot(data['pmra'][indices], data['pmdec'][indices], '.', markersize=3, color='red')
    plt.errorbar(data['pmra'], data['pmdec'], xerr=data['pmra_error'], yerr=data['pmdec_error'],
                 linestyle='', color='lightblue', marker='', zorder=0)
    plt.xlabel('pmra')
    plt.ylabel('pmdec')
    plt.xlim(-20, 20)
    plt.ylim(-20, 20)

    # (bp_rp,g)
    plt.subplot(223)
    # All stars by Gaia
    plt.plot(data['bp_rp'], data['g'], '.', markersize=1, color='gray')
    # Stars also in the Sampedro catalogue
    plt.plot(data['bp_rp'][indices], data['g'][indices], '.', markersize=3, color='red')
    plt.xlabel('BP-RP')
    plt.ylabel('G')
    plt.xlim(np.nanmin(data['bp_rp'][indices]), np.nanmax(data['bp_rp'][indices]))
    plt.ylim(np.nanmin(data['g'][indices]), np.nanmax(data['g'][indices]))

    # (parallax,g)
    plt.subplot(224)
    # All stars by Gaia
    plt.plot(data['parallax'], data['g'], '.', markersize=1, color='gray')
    # Stars also in the Sampedro catalogue
    plt.plot(data['parallax'][indices], data['g'][indices], '.', markersize=3, color='red')
    plt.errorbar(data['parallax'], data['g'], xerr=data['parallax_error'],
                 linestyle='', color='lightblue', marker='', zorder=0)
    plt.xlabel('Parallax')
    plt.ylabel('G')
    plt.xlim(0, np.nanmax(data['parallax'][indices]))
    plt.ylim(np.nanmin(data['g'][indices]), np.nanmax(data['g'][indices]))

    plt.savefig(data_folder + cluster_name + '_matched_plot.png')
    plt.close()


def plot_all_matches():
    try:
        with open('unverified_cluster_list.txt') as f:
            clusters = [c.strip() for c in f.readlines()]
    except FileNotFoundError:
        print("Gaia data not found. Make sure to retrieve it first.")

    i = 1
    total = len(clusters)
    for cluster in clusters:
        print(f"Plotting data of cluster {cluster}. Cluster {i}/{total}.")
        plot_match(cluster, 'verified.fits')
        i += 1


def main():
    #  get_unverified_cluster_data('sampedro_stars.fits', 'verified.fits')
    get_all_matches(2)
    plot_all_matches()


if __name__ == "__main__":
    main()
