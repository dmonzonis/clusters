import os
from astropy.io import fits
from astroquery.gaia import Gaia
import numpy as np
import matplotlib.pyplot as plt


dirname = os.path.dirname(__file__)


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


def get_unverified_clusters(identified_filename,
                            verified_filename,
                            dump_to_file=False):
    identified = fits.open(identified_filename)
    verified = fits.open(verified_filename)

    verified_set = {cluster[0] for cluster in verified[1].data}
    unverified_set = {cluster[0] for cluster in identified[1].data
                      if cluster[0] not in verified_set}

    identified.close()
    verified.close()

    return unverified_set


def main():
    unverified = get_unverified_clusters('sampedro_clusters.fit',
                                         'verified.fits',
                                         'unverified.fits')
    print(len(unverified))


if __name__ == "__main__":
    main()
