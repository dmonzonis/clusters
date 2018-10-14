from astropy.io import fits
from astroquery.gaia import Gaia

QUERY = "SELECT ra, dec, pmra, pmra_error, pmdec, pmdec_error, parallax, parallax_error, \
phot_g_mean_mag, bp_rp \
FROM gaiadr2.gaia_source \
WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),\
CIRCLE('ICRS',56.75,24.1167,2))=1 \
AND phot_g_mean_mag<17;"


def request_gaia_data(query=QUERY):
    job = Gaia.launch_job_async(QUERY, dump_to_file=True, verbose=True)
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
