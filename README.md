# Siluro-Devonian Communities
<a href="https://www.ucy.ac.cy/migrate/"><img align="right" src="www/MigratelogoShad.png" width="200" style="margin-top: 20px"></a>
Lithic and temporal variation predict composition of Siluro-Devonian vertebrate, invertebrate, and plant communities

<br>
<br>
<strong>lead investigator</strong>: Dr <a href="https://www.ucy.ac.cy/directory/en/profile/tmouts01">Theodora Moutsiou</a><br>
<strong>key personnel</strong>: Dr <a href="https://scholar.google.com.au/citations?user=BU25ogMAAAAJ&hl=en">Christian Reepmeyer</a>, Associate Professor <a href="https://www.ucy.ac.cy/directory/en/profile/demest">Stella Demesticha</a>, Prof <a href="https://www.presidency.gov.cy/cypresidency/cypresidency.nsf/All/FBA917CB5206BC95C225896B0023FAFD?OpenDocument">Vasiliki Kassianidou</a>, Dr <a href="https://www.cut.ac.cy/faculties/fet/ceg/staff/athos.agapiou/?languageId=1">Athos Agapiou</a>, Dr <a href="https://www.researchgate.net/profile/Zomenia-Zomeni">Zomenia Zomeni</a>, Professor <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey Bradshaw</a><br>
<strong>collaborators</strong>: Dr <a href="https://globalecologyflinders.com/people/#COORDINATOR">Frédérik Saltré</a>, Dr <a href="https://qcnr.usu.edu/directory/envs/faculty/crabtree-stefani">Stefani Crabtree</a>
<br>

## Abstract
While vertebrate communities of the Siluro-Devonian have been well-described, co-variation in invertebrate and plant communities from the same period has not been adequately assessed. Using methods from community ecology, we tested how the combined vertebrate, invertebrate, and plant communities from 105 fossiliferous sites around the world varied relative to sedimentology and time. We built logistic linear models to examine variation in the three main communities (vertebrates, invertebrates, plants). We found that lithic type (6 categories) explained the most variation (% deviance explained [%DE] = 8.0), followed by geological period (4 categories); however, community type had little explanatory power (%DE = 0.05). Lithic variation in Silurian communities was driven predominately by low presence probability in limestone, whereas variation in Devonian communities was driven mainly by higher presence probabilities in sandstone and shale. There was little support for a period×lithics interaction, nor was much variation explained by any factor when we divided communities by major palaeo-environmental category (marine, alluvial/deltaic, freshwater/estuarine). We then applied a stochastic variant of a permutation analysis of variance to the entire community (10 vertebrate, 20 invertebrate, 7 plant taxa) to examine the relative explanatory power of lithic and temporal variation on composition. There was an effect of lithics (R2perm = 0.166; pperm = 0.0094; driven mainly by limestone) and categorical period (R2perm = 0.271; pperm = 0.0179; driven mainly by the Middle Devonian) (but not their interaction: pperm = 0.152) on community composition. Treating time as a continuous variable (age) also demonstrated an effect on community composition (R2perm = 0.093; pperm = 0.0052), with incidence tending to rise in most taxa from 435 to 360 million years ago. Model predictions can now be used to show spatial variation in the relative presence probability for any of the taxa we modelled provided detailed maps of lithic distribution are available.

## <a href="https://github.com/cjabradshaw/SiluroDevonianCommunities/tree/main/scripts">Scripts</a>
- <code>silurodev.R</code>: main code to generate results
- <code>r.squared.R</code>: functions to estimate goodness of fit for linear models
- <code>new_lmer_AIC_tables3.r</code>: functions to compare linear models

## <a href="https://github.com/cjabradshaw/SiluroDevonianCommunities/tree/main/data">Data</a>
- <em>sildevcomp2.txt</em>: fossiliferous site database

## R (v4.3.2) libraries
<code>performance</code>, <code>sjPlot</code>, <code>lme4</code>, <code>ggplot2</code>, <code>stringr</code>, <code>vegan</code>, <code>parallel</code>, <code>vcdExtra</code>, <code>plyr</code>, <code>dplyr</code>

<br>
<p><a href="https://www.ucy.ac.cy"><img align="bottom-left" src="www/UCypruslogo.png" alt="UCyprus logo" height="40" style="margin-top: 20px"></a> &nbsp; <a href="http://www.dainst.org"><img align="bottom-left" src="www/DAIlogo.png" alt="DAI logo" height="55" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" height="30" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp-2.png" alt="GEL logo" height="55" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://EpicAustralia.org.au"><img align="bottom-left" src="www/CabahFCL.jpg" alt="CABAH logo" height="40" style="margin-top: 20px"></a> &nbsp; <a href="https://www.cut.ac.cy"><img align="bottom-left" src="www/CUTlogoblack.png" alt="CUT logo" height="50" style="margin-top: 20px"></a><a href="https://www.moa.gov.cy/moa/gsd/gsd.nsf/dmlIndex_en/dmlIndex_en"><img align="bottom-left" src="www/CGSlogo.png" alt="CGS logo" height="45" style="margin-top: 20px"></a></p>

