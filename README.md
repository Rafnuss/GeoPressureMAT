# GeoPressureMAT

Matlab code to analyse geolocator using atmospheric pressure.

> **Warning**
> This repository will not be maintain. Check out the R package [GeoPressureR](https://raphaelnussbaumer.com/GeoPressureR/).

This repository accompagny the studies:

> Nussbaumer, R., Gravey, M., Briedis, M., & Liechti, F. (2022). Global positioning with animal-borne pressure sensors. Methods in Ecology and Evolution, 00, 1– 14. <https://doi.org/10.1111/2041-210X.14043>

> Nussbaumer, R., Gravey, M., Briedis, M., Liechti, F. & Sheldon, D. (2022) Reconstructing bird trajectories from pressure and wind data using a highly optimised hidden Markov model. PREPRINT (Version 2). <https://doi.org/10.21203/rs.3.rs-1693751/v2>

The [GeoPressureAPI](https://github.com/Rafnuss/GeoPressureAPI) is a JSON API that makes it easy to compute the mismatch of a geolocator pressure timeserie with the atmospheric pressure from ERA5-LAND reanalysis data. This can be used with any computational language.

The [Google Earth Engine App](https://rafnuss.users.earthengine.app/view/pressuregeolocator) visualize the position at each the stationary periods with covariance such as landcover and NDVI.
