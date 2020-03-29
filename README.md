# Atmospheric convection plays a key role in the climate of tidally-locked terrestrial exoplanets: insights from high-resolution simulations
**Denis E. Sergeev, F. Hugo Lambert, Nathan J. Mayne, Ian Boutle, James Manners, and Krisztian Kohary**

**2020**

**Submitted to *Astrophysical Journal***

Code used to process and visualise the model output.
Model data are available upon request (raw data O(100 Gb)).

Notebooks for each individual figure as well as for two data tables are in the [`code/` directory](code), while the figures themselves are in the `plots/` directory.

| Figure                                                                                                         | Notebook                                                                              | Dependencies                                  |
|:---------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------|:----------------------------------------------|
| Overview of the model set-up for the Trappist-1e case (interactive 3D Figure)                                  | [Fig01-Model-Set-Up-3D.ipynb](code/Fig01-Model-Set-Up-3D.ipynb)                       | [pp_ns_data.py](code/pp_ns_data.py)           |
| Surface temperature and horizontal wind vectors in the upper troposphere                                       | [Fig02-Map-Sfc-Temp-Winds.ipynb](code/Fig02-Map-Sfc-Temp-Winds.ipynb)                 | [pp_main_gl_data.py](code/pp_main_gl_data.py) |
| Time average vertical profiles of temperature and water vapor at the sub-stellar point and its antipode        | [Fig03-Vert-Prof-Temp-Hum.ipynb](code/Fig03-Vert-Prof-Temp-Hum.ipynb)                 | [pp_main_gl_data.py](code/pp_main_gl_data.py) |
| Eddy geopotential height and eddy components of horizontal wind vectors 250 hPa and air temperature at 700 hPa | [Fig04-Map-Gpt-Hgt-Temp.ipynb](code/Fig04-Map-Gpt-Hgt-Temp.ipynb)                     | [pp_main_gl_data.py](code/pp_main_gl_data.py) |
| Atmospheric circulation regimes                                                                                | [Fig05-Rossby-Rhines.ipynb](code/Fig05-Rossby-Rhines.ipynb)                           | [pp_main_gl_data.py](code/pp_main_gl_data.py) |
| Distribution of clouds and precipitation for the Trappist-1e and Proxima b cases                               | [Fig06-Map-Cloud-Precip.ipynb](code/Fig06-Map-Cloud-Precip.ipynb)                     | [pp_main_gl_data.py](code/pp_main_gl_data.py) |
| Meridional average of vertically integrated divergence of MSE and its components                               | [Fig07-MSE-Flux-Divergence.ipynb](code/Fig07-MSE-Flux-Divergence.ipynb)               | [pp_main_gl_data.py](code/pp_main_gl_data.py) |
| Latitudinal cross-section of the zonal transport of sensible and latent heat across the eastern terminator     | [Fig08-Vert-Cross-Heat-Flux.ipynb](code/Fig08-Vert-Cross-Heat-Flux.ipynb)             | [pp_main_gl_data.py](code/pp_main_gl_data.py) |
| Snapshot of top-of-atmosphere outgoing longwave radiation in the Trappist-1e simulation                        | [Fig09-Map-TOA-OLR-Trap1e.ipynb](code/Fig09-Map-TOA-OLR-Trap1e.ipynb)                 |                                               |
| Snapshot of top-of-atmosphere outgoing longwave radiation in the Proxima b simulation                          | [Fig10-Map-TOA-OLR-Proxb.ipynb](code/Fig10-Map-TOA-OLR-Proxb.ipynb)                   |                                               |
| Histograms of instantaneous TOA OLR values in the substellar region                                            | [Fig11-Hist-TOA-OLR.ipynb](code/Fig11-Hist-TOA-OLR.ipynb)                             |                                               |
| Snapshot of convection in the substellar region in the Proxima b case at the end of the HighRes simulation     | [Fig12-Map-Precip-Vert-Wind.ipynb](code/Fig12-Map-Precip-Vert-Wind.ipynb)             |                                               |
| Vertical profiles of the cloud liquid water and ice mixing ratio in the substellar region                      | [Fig13-Vert-Prof-Cloud-Condensate.ipynb](code/Fig13-Vert-Prof-Cloud-Condensate.ipynb) | [pp_ns_data.py](code/pp_ns_data.py)           |
|                                                                                                                |                                                                                       | [ns_mean_vprof.py](code/ns_mean_vprof.py)     |
| Vertical profiles of temperature increments due to the latent heating in the substellar region                 | [Fig14-Vert-Prof-Latent-Heating.ipynb](code/Fig14-Vert-Prof-Latent-Heating.ipynb)     | [pp_ns_data.py](code/pp_ns_data.py)           |
|                                                                                                                |                                                                                       | [ns_mean_vprof.py](code/ns_mean_vprof.py)     |
| Day-night surface temperature difference vs wind divergence in the free troposphere of the substellar region   | [Fig15-Day-Night-Impact-Lin-Reg.ipynb](code/Fig15-Day-Night-Impact-Lin-Reg.ipynb)     | [pp_ens_gl_data.py](code/pp_ens_gl_data.py)   |
|                                                                                                                |                                                                                       | [aggr_ens_output.py](code/aggr_ens_output.py) |
|                                                                                                                |                                                                                       | [pp_ns_data.py](code/pp_ns_data.py)           |
|                                                                                                                |                                                                                       | [aggr_ns_data.py](code/aggr_ns_data.py)       |
| Mean global, day-side, and night-side surface temperature in the control and sensitivity simulations           | [Tab03-Mean-Sfc-Temp.ipynb](code/Tab03-Mean-Sfc-Temp.ipynb)                           | [pp_main_gl_data.py](code/pp_main_gl_data.py) |
| Global top-of-atmosphere cloud radiative effect in the control and sensitivity simulations                     | [Tab04-CRE.ipynb](code/Tab04-CRE.ipynb)                                               | [pp_main_gl_data.py](code/pp_main_gl_data.py) |

`page/` directory is for hosting an HTML page with an [interactive version of Fig. 1](https://github.com/dennissergeev/exoconvection-apj-2020).
