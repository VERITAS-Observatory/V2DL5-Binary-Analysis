"""
Target definition.
"""

import logging

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord, name_resolve
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

        self.exclusion_mask = self.get_exclusion_mask(
            on_region_dict=args_dict["datasets"]["on_region"],
            on_region_exclusion_radius=args_dict["datasets"]["exclusion_region"]["radius"],
        )

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

    def get_exclusion_mask(self, on_region_dict=None, on_region_exclusion_radius=None):
        """
        Defines a mask for the exclusion regions.

        Parameters
        ----------
        on_region : dict
            on_region dictionary.

        """

        exclusion_regions = []

        # on region
        if on_region_dict is not None and on_region_exclusion_radius is not None:
            exclusion_regions.append(
                CircleSkyRegion(
                    center=self.get_target(sky_coord=on_region_dict, info=False),
                    radius=Angle(on_region_exclusion_radius),
                )
            )
            self._logger.info(f"On region exclusion: {exclusion_regions[-1]}")

        # bright stars
        # TODO

        # exclusion mask
        geom = WcsGeom.create(
            npix=(150, 150), binsz=0.05, skydir=self.target.galactic, proj="TAN", frame="icrs"
        )

        self._logger.info("Number of exclusion regions: %d", len(exclusion_regions))
        return ~geom.region_mask(exclusion_regions)
