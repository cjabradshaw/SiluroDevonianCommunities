# Siluro-Devonian Communities
<img align="right" src="www/ageyr perm.jpg" width="400" style="margin-top: 20px"></a>
Lithofacies and temporal variation predict composition of Siluro-Devonian vertebrate, invertebrate, and plant communities
<br>
<br>
lead: Professor <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey Bradshaw</a><br>
collaborators: <a href="https://www.uqar.ca/universite/a-propos-de-l-uqar/departements/departement-de-biologie-chimie-et-geographie/cloutier-richard">Richard Cloutier</a>, <a href="https://prc.msu.ac.th/eng/personnel-palaeontological-research-and-education-centre-msu/">Clive Burrett</a>, <a href="https://prc.msu.ac.th/eng/personnel-palaeontological-research-and-education-centre-msu/">Mongkol Udchachon</a>
<br>
## Abstract
While Siluro-Devonian vertebrate communities have been well-described, co-variation in invertebrate and plant communities from the same period has not been adequately assessed. Using methods from community ecology, we tested how the combined vertebrate, invertebrate, and plant communities from 105 fossiliferous sites around the world varied relative to lithofacies and time. We built logistic linear models to examine variation in the three main communities (vertebrates, invertebrates, plants). We found that the type of lithofacies (6 categories) explained the most variation (% deviance explained [%DE] = 8.0), followed by geological period (4 categories); however, community type had little explanatory power (%DE = 0.05). Variation among lithofacies in Silurian communities was driven predominately by low presence probability in limestone, whereas variation in Devonian communities was driven mainly by higher presence probabilities in sandstone and shale. There was little support for a period×lithofacies interaction, nor was much variation explained by any factor when we divided communities by major paleo-environmental category (marine, alluvial/deltaic, freshwater/estuarine). We then applied a stochastic variant of a permutation analysis of variance to the entire community (10 vertebrate, 20 invertebrate, 7 plant taxa) to examine the relative explanatory power of lithofacies and temporal variation on composition. There was an effect of the lithofacies (R<sup>2</sup><sub>perm</sub> = 0.166; <em>p</em><sub>perm</sub>; driven mainly by limestone) and categorical period (R<sup>2</sup><sub>perm</sub> = 0.271; <em>p</em><sub>perm</sub> = 0.0179; driven mainly by the Middle Devonian) (but not their interaction: <em>p</em><sub>perm</sub> = 0.152) on community composition. Treating time as a continuous variable (age) also demonstrated an effect on community composition (R<sup>2</sup><sub>perm</sub> = 0.093; <em>p</em><sub>perm</sub> = 0.0052), with incidence tending to rise in most taxa from 435 to 360 million years ago. Model predictions can now be used to show spatial variation in the relative presence probability for any of the taxa we modelled provided detailed distribution maps of lithofacies are available.

## <a href="https://github.com/cjabradshaw/SiluroDevonianCommunities/tree/main/scripts">Scripts</a>
- <code>silurodev.R</code>: main code to generate results
- <code>r.squared.R</code>: functions to estimate goodness of fit for linear models
- <code>new_lmer_AIC_tables3.r</code>: functions to compare linear models
- <code>palaeocoords.ipynb</code>: Jupyter Notebook with Python code to project current-era lat/lon coordinates to palaeo-coordinates, with corresponding plate tectonics, for anytime in the last 1 billion years (modified code supplied by @<a href="https://github.com/amer7632">amer7632</a>; detailed instructions <a href="https://github.com/GPlates/gplately/blob/master/Notebooks/03-WorkingWithPoints.ipynb">here</a>)

## <a href="https://github.com/cjabradshaw/SiluroDevonianCommunities/tree/main/data">Data</a>
- <em>sildevcomp2.txt</em>: fossiliferous site database
- <a href="https://github.com/cjabradshaw/SiluroDevonianCommunities/tree/main/data/plate-model-repo/">plate-model-repo</a>: data for hindcasted tectonic reconstruction and generation of palaeo coordinates (maintain sub-folder hierarchy & unzip files in 'Topologies' sub-folder)
- <a href="https://github.com/cjabradshaw/SiluroDevonianCommunities/tree/main/data/sitecoords/">sitecoords</a>: .csv files with site coordinates — all: <em>sitelocs.csv</em>; and divided by major period: Silurian (SIL), Lochkovian (LOCH), Emsian (EMS), Eifelian (EIF), Givetian (GIV), Frasnian (FRAS), Famennian (FAM)
<img align="right" src="www/FAMsites.jpg" width="400" style="margin-top: 20px"></a>

## R libraries
<code>performance</code>, <code>sjPlot</code>, <code>lme4</code>, <code>ggplot2</code>, <code>stringr</code>, <code>vegan</code>, <code>parallel</code>, <code>vcdExtra</code>, <code>plyr</code>, <code>dplyr</code>

## GPlately
Install <code><a href="https://github.com/GPlates/gplately?tab=readme-ov-file">GPlately</a></code>, a Python interface to accelerate spatio-temporal data analysis leveraging <code><a href="https://www.gplates.org/docs/pygplates/index.html">pyGPlates</a></code> and <code><a href="https://github.com/EarthByte/PlateTectonicTools">PlateTectonicTools</a></code>
(references: <a href="https://doi.org/10.1016/j.earscirev.2020.103477">Merdith et al.</a> 2021 <em>Earth-Science Reviews</em> 214:103477; <a href="https://doi.org/10.1002/gdj3.185">Mather et al.</a> 2024 <em>Geosciences Data Journal</em> 11:3-10)

## Python dependencies
<code><a href="https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation">pyGPlates</a></code>, <code><a href="https://pypi.org/project/plate-model-manager/">plate-model-manager</a></code>, <code><a href="https://shapely.readthedocs.io/en/stable/project.html#installing-shapely">Shapely</a></code>, <code><a href="https://numpy.org/install/">NumPy</a></code>. <code><a href="https://scipy.org/install/">SciPy</a></code>, <code><a href="https://matplotlib.org/stable/users/installing/index.html">Matplotlib</a></code>, <code><a href="https://scitools.org.uk/cartopy/docs/latest/index.html#getting-started">Cartopy</a></code>, <code><a href="https://github.com/fatiando/pooch">Pooch</a></code>, <code><a href="https://geopandas.org/en/stable/getting_started.html">GeoPandas</a></code>, <code><a href="https://unidata.github.io/netcdf4-python/#quick-install">netCDF4</a></code>

<br>
<a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" height="40" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp-2.png" alt="GEL logo" height="85" style="margin-top: 20px"></a> &nbsp; <a href="https://www.uqar.ca"><img align="bottom-left" src="www/UQARlogo.png" alt="CUT logo" height="70" style="margin-top: 20px"></a><a href="https://www.msu.ac.th/eng/"><img align="bottom-left" src="www/MahasarakhamUlogo.png" alt="CGS logo" height="95" style="margin-top: 20px"></a></p>

