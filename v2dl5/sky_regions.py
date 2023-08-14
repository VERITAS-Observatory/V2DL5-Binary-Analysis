"""
Sky regions definition.
"""

import logging

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord, name_resolve
from astropy.io import fits
from astropy.table import Table
from gammapy.maps import WcsGeom
from regions import CircleSkyRegion


class SkyRegions:
    """
    Define sky regions required for the analysis:
    - target coordinates
    - on region
    - exclusion regions

    Parameters
    ----------
    args_dict : dict
        Dictionary of configuration arguments.

    """

    def __init__(self, args_dict=None):
        self._logger = logging.getLogger(__name__)

        self.target = self.get_target(sky_coord=args_dict["observations"]["obs_cone"])
        self.on_region = self.get_on_region(on_region_dict=args_dict["datasets"]["on_region"])
        self.exclusion_mask = self.get_exclusion_mask(args_dict)

    def get_target(self, sky_coord=None, info=True):
        """
        Defines a SkyCoord object for the target.
        Reads target coordinates from Simbad if target name is given.

        Parameters
        ----------
        sky_coord : dict
            sky_coord dictionary.

        info: bool
            Print info about target.

        """

        target = None
        if sky_coord.get("target", None) is None:
            if sky_coord.get("frame", None) == "icrs":
                target = SkyCoord(
                    ra=sky_coord.get("lon", None),
                    dec=sky_coord.get("lat", None),
                    frame=sky_coord.get("frame", None),
                )
            elif sky_coord.get("frame", None) == "galactic":
                target = SkyCoord(
                    l=sky_coord.get("lon", None),
                    b=sky_coord.get("lat", None),
                    frame=sky_coord.get("frame", None),
                )
            else:
                raise ValueError("Unsupported coordinate frame")
            self._logger.info(f"Target coordinates given: {target}")
        else:
            try:
                target = SkyCoord.from_name(sky_coord["target"])
                self._logger.info("Target %s found in Simbad.", sky_coord["target"])
            except name_resolve.NameResolveError:
                self._logger.error('Target "%s" not found in Simbad.', sky_coord["target"])
                raise

        if info:
            self._logger.info("Target name: %s", sky_coord.get("target", "target name not given"))
            self._logger.info("Target coordinates: %s", target)

        return target

    def get_on_region(self, on_region_dict=None):
        """
        Defines a CircleSkyRegion object for the on region.

        Parameters
        ----------
        on_region : dict
            on_region dictionary.

        Returns
        -------
        on_region : CircleSkyRegion
            on region.

        """

        on_region = CircleSkyRegion(
            center=self.get_target(sky_coord=on_region_dict, info=False),
            radius=Angle(on_region_dict.get("radius", 0.5 * u.deg)),
        )

        self._logger.info(f"On region: {on_region}")

        return on_region

    def get_exclusion_mask(self, args_dict):
        """
        Defines a mask for the exclusion regions.

        Parameters
        ----------
        args_dict: dict
            Dictionary of configuration arguments.

        """

        exclusion_regions = []

        # on region
        on_region_dict = args_dict["datasets"]["on_region"]
        on_region_exclusion_radius = args_dict["datasets"]["exclusion_region"]["on_radius"]
        if on_region_dict is not None and on_region_exclusion_radius is not None:
            exclusion_regions.append(
                CircleSkyRegion(
                    center=self.get_target(sky_coord=on_region_dict, info=False),
                    radius=Angle(on_region_exclusion_radius),
                )
            )
            self._logger.info(f"On region exclusion: {exclusion_regions[-1]}")

        # bright star exclusion
        exclusion_regions.extend(
            self._read_bright_star_catalogue(
                exclusion_region_dict=args_dict["datasets"]["exclusion_region"]
            )
        )

        # exclusion mask
        geom = WcsGeom.create(
            npix=(150, 150), binsz=0.05, skydir=self.target.galactic, proj="TAN", frame="icrs"
        )

        self._logger.info("Number of exclusion regions: %d", len(exclusion_regions))
        return ~geom.region_mask(exclusion_regions)

    def _read_bright_star_catalogue(self, exclusion_region_dict=None):
        """
        Read bright star catalogue from file.

        Parameters
        ----------
        exclusion_region_dict: dict
            Dictionary of configuration arguments.

        Returns
        -------
        List of CircleSkyRegion
            List exclusion regions due to bright stars.

        """

        _exclusion_regions = []
        if exclusion_region_dict["star_file"] is None:
            return _exclusion_regions

        self._logger.info(
            "Reading bright star catalogue from %s", exclusion_region_dict["star_file"]
        )
        hip = fits.open(exclusion_region_dict["star_file"])
        catalogue = Table(hip[1].data)
        catalogue = catalogue[
            catalogue["Vmag"] + catalogue["B-V"] < exclusion_region_dict["magnitude_B"]
        ]

        angular_separation = self.target.separation(
            SkyCoord(catalogue["_RA_icrs"], catalogue["_DE_icrs"], unit=u.deg)
        )
        catalogue["angular_separation"] = angular_separation

        catalogue = catalogue[
            catalogue["angular_separation"] < u.Quantity(exclusion_region_dict["fov"])
        ]

        self._logger.info(
            f"Number of stars in the catalogue passing cuts on magnitude and FOV: {len(catalogue)}"
        )
        print(catalogue.pprint_all())

        for row in catalogue:
            _exclusion_regions.append(
                CircleSkyRegion(
                    center=SkyCoord(row["_RA_icrs"], row["_DE_icrs"], unit=u.deg),
                    radius=Angle(exclusion_region_dict["star_exclusion_radius"]),
                )
            )

        return _exclusion_regions
