# Siluro-Devonian Communities
<a href="https://www.ucy.ac.cy/migrate/"><img align="right" src="www/MigratelogoShad.png" width="200" style="margin-top: 20px"></a>
Lithic and temporal variation predict composition of Siluro-Devonian vertebrate, invertebrate, and plant communities
<br>
<br>
lead: Professor <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey Bradshaw</a><br>
collaborators: <a href="https://www.uqar.ca/universite/a-propos-de-l-uqar/departements/departement-de-biologie-chimie-et-geographie/cloutier-richard">Richard Cloutier</a>, <a href="https://prc.msu.ac.th/eng/personnel-palaeontological-research-and-education-centre-msu/">Clive Burrett</a>, <a href="https://prc.msu.ac.th/eng/personnel-palaeontological-research-and-education-centre-msu/">Mongkol Udchachon</a>
<br>
## Abstract
While vertebrate communities of the Siluro-Devonian have been well-described, co-variation in invertebrate and plant communities from the same period has not been adequately assessed. Using methods from community ecology, we tested how the combined vertebrate, invertebrate, and plant communities from 105 fossiliferous sites around the world varied relative to sedimentology and time. We built logistic linear models to examine variation in the three main communities (vertebrates, invertebrates, plants). We found that lithic type (6 categories) explained the most variation (% deviance explained [%DE] = 8.0), followed by geological period (4 categories); however, community type had little explanatory power (%DE = 0.05). Lithic variation in Silurian communities was driven predominately by low presence probability in limestone, whereas variation in Devonian communities was driven mainly by higher presence probabilities in sandstone and shale. There was little support for a period×lithics interaction, nor was much variation explained by any factor when we divided communities by major palaeo-environmental category (marine, alluvial/deltaic, freshwater/estuarine). We then applied a stochastic variant of a permutation analysis of variance to the entire community (10 vertebrate, 20 invertebrate, 7 plant taxa) to examine the relative explanatory power of lithic and temporal variation on composition. There was an effect of lithics (R2perm = 0.166; pperm = 0.0094; driven mainly by limestone) and categorical period (R2perm = 0.271; pperm = 0.0179; driven mainly by the Middle Devonian) (but not their interaction: pperm = 0.152) on community composition. Treating time as a continuous variable (age) also demonstrated an effect on community composition (R2perm = 0.093; pperm = 0.0052), with incidence tending to rise in most taxa from 435 to 360 million years ago. Model predictions can now be used to show spatial variation in the relative presence probability for any of the taxa we modelled provided detailed maps of lithic distribution are available.

## <a href="https://github.com/cjabradshaw/SiluroDevonianCommunities/tree/main/scripts">Scripts</a>
- <code>silurodev.R</code>: main code to generate results
- <code>r.squared.R</code>: functions to estimate goodness of fit for linear models
- <code>new_lmer_AIC_tables3.r</code>: functions to compare linear models

## <a href="https://github.com/cjabradshaw/SiluroDevonianCommunities/tree/main/data">Data</a>
- <em>sildevcomp2.txt</em>: fossiliferous site database

## R libraries
<code>performance</code>, <code>sjPlot</code>, <code>lme4</code>, <code>ggplot2</code>, <code>stringr</code>, <code>vegan</code>, <code>parallel</code>, <code>vcdExtra</code>, <code>plyr</code>, <code>dplyr</code>

<br>
<a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" height="40" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp-2.png" alt="GEL logo" height="65" style="margin-top: 20px"></a> &nbsp; <a href="https://www.uqar.ca"><img align="bottom-left" src="www/UQARlogo.png" alt="CUT logo" height="50" style="margin-top: 20px"></a><a href="https://www.msu.ac.th/eng/"><img align="bottom-left" src="www/MahasarakhamUlogo.png" alt="CGS logo" height="55" style="margin-top: 20px"></a></p>

