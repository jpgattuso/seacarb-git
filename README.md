seacarb
=======

Carbonate chemistry with R

In 2003, Aurélien Proye and I put together seacarb, an R package that calculates various parameters of the carbonate system in seawater. The package was subsequently upgraded. In 2008, a new version (2.0) was built with the assistance of Héloïse Lavigne. Five additional functions were included to assist the design of perturbation experiments to investigate ocean acidification and other functions were revised in order to strictly follow the "Guide to best practices for ocean CO2 measurements" (Dickson et al., 2007). Seacarb uses equations mostly from the following sources:

- Dickson A. G., Sabine C. L. & Christian J. R., 2007. Guide to best practices for ocean CO2 measurements. PICES Special Publication 3:1-191.
- DOE, 1994. Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water; version 2. Dickson, A. G. and Goyet, C., editors. ORNL/CDIAC-74, 1994.
- Frankignoulle M., 1994. A complete set of buffer factors for acid/base CO2 system in seawater. Journal of Marine Systems 5: 111-118.
- Orr J. C., Epitalon J.-M., Dickson A. G. & Gattuso J.-P., 2018. Routine uncertainty propagation for the marine carbon dioxide system. Marine Chemistry 207:84-107.
- Zeebe R. E. & Wolf-Gladrow D. A., 2001. CO2 in seawater: equilibrium, kinetics, isotopes. Amsterdam: Elsevier, 346 pp.

A portion of the code has been adapted, with permission from the authors, from the Matlab files mentioned above. I am grateful to Richard Zeebe and Dieter Wolf-Gladrow for that. Portions of code and/or corrections have also been contributed by Jean-Marie Epitalon (2004, 2018), Bernard Gentili (2006), Karline Soetaert (2007) and Jim Orr (2007, 2010, 2018). Héloïse Lavigne contributed a lot to seacarb in 2008 and 2009; this was critical to the launch of version 2.0 and subsequent updates. Jean-Marie Epitalon considerably improved the code, leading to much faster calculations in version 3.0.

This program is provided free under the GNU General Public License (GNU GPL). The current version (3.0) will be improved using the comments that I will receive. If you are new to R, please check the manuals and FAQs available on the R-project web site to get information on how to install R and the seacarb package on your system. Please only report and comment on seacarb, not on general problems related to R.

Briefly, after installing R and if you have an Internet connection, here is the simplest way to install seacarb:

- Launch R
- To install seacarb (to be done only once), type the following command: install.packages("seacarb")
- To load the seacarb package into memory in order to use it (to be done each time R is launched), type the following command: library(seacarb)

The seacarb package can be downloaded from the Comprehensive R Archive Network (CRAN). The documentation is included in the package and is accessible using standard R commands. Please give due credit to the publications mentioned above and cite seacarb as follows:

Gattuso J.-P., Epitalon J.-M., Lavigne H. & Orr J., 2018. seacarb: seawater carbonate chemistry. R package version 3.2.10. http://CRAN.R-project.org/package=seacarb

