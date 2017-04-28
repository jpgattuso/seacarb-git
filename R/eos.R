# Copyright (C) 2008 Jean-Marie Epitalon and Heloise Lavigne and Aurelien Proye and Jean-Pierre Gattuso  
# with a most valuable contribution of Bernard Gentili
#
# This file is part of seacarb.
#/
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
sp2sa_chem <- function (SP, TA=2300e-6, DIC=2000e-6, NO3=0, SIOH4=0)
#
# This function converts from Practical (SP) to Absolute Salinity (SA)
# based on Total Alkalinity, Dissolved Inorganic Carbon 
# and Nitrate and Silicate concentrations.
#
# It is made of two parts :
#   1) conversion from Practical to Reference Salinity (major part)
#   2) conversion from Reference to Absolute Salinity (Absolute Anomaly, minor part)
#
{
    # DIC and TA values of standard sea water
    TA_stdsw  <- 0.0023  * SP /35
    DIC_stdsw <- 0.00208 * SP /35
    
    # Absolute salinity anomaly (g/kg)
    sal_anomaly = 55.6 * (TA - TA_stdsw) + 4.7 * (DIC - DIC_stdsw) + 38.9 * NO3 + 50.7 * SIOH4
    
    SR <- SP * 35.16504 / 35
    SA <- SR + sal_anomaly
    return(SA)
}


sp2sa_geo <- function (SP, P=0, long=1.e20, lat=1.e20)
#
# This function converts from Practical (SP) to Absolute Salinity (SA)
# based on Pressure P in dbar and geographic location : 
# latitude in decimal degrees [-90 ... 90] 
# and longitude in decimal degrees [ 0 ... +360 ] or  [ -180 ... +180 ]
#
# Latitude and longitude are optional. 
# If any is omitted, an arbitrary geographic point is chosen : central equatorial Atlantic.
# Note that this implies an error on computed SA ranging from 0 up to 0.02 g/kg
{
    if ((length(lat) == 1 && lat == 1.e20) || (length(long) == 1 && long == 1.e20))
    {
        lat <- 0; long <- -25;
    }
    return (gsw_SA_from_SP(SP,P,long,lat))
}


sa2sp_chem <- function (SA, TA=2300e-6, DIC=2000e-6, NO3=0, SIOH4=0)
#
# This function does the reverse of above function SA_from_SP_chem()
# It converts from Absolute (SA) to Practical Salinity (SP)
# based on chemical composition of sea water : Total Alkalinity, Dissolved Inorganic Carbon 
# and Nitrate and Silicate concentrations.
#
# Reverse conversion (SA from SP) follows this equation :
#  SP = 1/r * (SA - sa_anomaly)
#    with
#        r = 35.16504 / 35
#  SP = 1/r * (SA - 55.6 * (TA - x*SP) - 4.7 * (DIC - y*SP) - nutrients)
#    with
#       x  = 0.00230 / 35
#       y  = 0.00208 / 35
#       nutrients = 38.9 * NO3 + 50.7 * SIOH4
#
#  SP * (1 - 55.6 * x/r - 4.7 * y/r) = 1/r * (SA - 55.6 * TA - 4.7 * DIC - nutrients)
#
#  SP * (1 - 55.6 * 0.00230 / 35.16504 - 4.7 * 0.00208 / 35.16504) = ....
#
#  SP * 0.9960854303 = 1/r * (SA - 55.6 * TA - 4.7 * DIC - nutrients)
#
#  SP = 0.999218211671 * (SA - 55.6 * TA - 4.7 * DIC - nutrients)
{
    SP <- 0.999218211671 * (SA - 55.6 * TA - 4.7 * DIC - 38.9 * NO3 - 50.7 * SIOH4)
    return(SP)
}


sa2sp_geo <- function (SA, P=0, long=1.e20, lat=1.e20)
#
# This function converts from Absolute (SA) to Practical Salinity (SP)
# based on Pressure P in dbar and geographic location : 
#   - latitude in decimal degrees [-90 ... 90] 
#   - longitude in decimal degrees [ 0 ... +360 ] or  [ -180 ... +180 ]
#
# Latitude and longitude are optional. 
# If any is omitted, an arbitrary geographic point is chosen : central equatorial Atlantic.
# Note that this implies an error on computed SP ranging from 0 up to 0.02 g/kg
{
    if ((length(lat) == 1 && lat == 1.e20) || (length(long) == 1 && long == 1.e20))
    {
        lat <- 0; long <- -25;
    }
    return (gsw_SP_from_SA(SA,P,long,lat))
}


eos2teos_geo <- function (SP, T, P=0, long=1.e20, lat=1.e20)
#
# This function converts Temperature and Salinity from the EOS-80 to the TEOS-10 standard.
# 
# Input parameters :
#   - SP       :   Practical salinity in psu
#   - T        :   In situ Temperature in degrees C
#   - P        :   Sea Pressure in dbar
#   - long     :   longitude in decimal degrees [ 0 ... +360 ] or  [ -180 ... +180 ]
#   - lat      :   latitude in decimal degrees [-90 ... 90] 
#
# Latitude and longitude are optional. 
# If any is omitted, an arbitrary geographic point is chosen : central equatorial Atlantic.
# Note that this implies an error on computed SP ranging from 0 up to 0.02 g/kg
#
# Output : a data frame with the following colums
#    - CT      :   Conservative Temperature
#    - SA      :   Absolute Salinity
#
{
    if ((length(lat) == 1 && lat == 1.e20) || (length(long) == 1 && long == 1.e20))
    {
        lat <- 0; long <- -25;
    }
    # convert salinity
    SA <- gsw_SA_from_SP(SP, P, long, lat)
    # convert temperature
    CT <- gsw_CT_from_t (SA, T, P)
    
    return (as.data.frame(cbind(CT, SA)))
}

teos2eos_geo <- function (SA, CT, P=0, long=1.e20, lat=1.e20)
#
# This function converts Temperature and Salinity from the TEOS-10 to the EOS-80 standard.
# 
# Input parameters :
#   - SA       :   Absolute salinity in g/kg
#   - CT       :   Conservative Temperature in degrees C
#   - P        :   Sea Pressure in dbar
#   - long     :   longitude in decimal degrees [ 0 ... +360 ] or  [ -180 ... +180 ]
#   - lat      :   latitude in decimal degrees [-90 ... 90] 
#
# Latitude and longitude are optional. 
# If any is omitted, an arbitrary geographic point is chosen : central equatorial Atlantic.
# Note that this implies an error on computed SP ranging from 0 up to 0.02 g/kg
#
# Output : a data frame with the following colums
#    - T       :   in situ Temperature
#    - SP      :   Practical Salinity
#
{
    if ((length(lat) == 1 && lat == 1.e20) || (length(long) == 1 && long == 1.e20))
    {
        lat <- 0; long <- -25;
    }
    # convert temperature
    T <- gsw_t_from_CT (SA, CT, P)
    # convert salinity
    SP <- gsw_SP_from_SA(SA, P, long, lat)
    
    return (as.data.frame(cbind(T, SP)))
}


eos2teos_chem <- function (SP, T, P=0, TA=2300e-6, DIC=2000e-6, NO3=0, SIOH4=0)
#
# This function from Temperature and Salinity from EOS-80 to TEOS-10 standard.
# 
# Input parameters :
#   - SP       :   Practical salinity in psu
#   - T        :   in situ Temperature in degrees C
#   - TA       :   Total Alkalinity in mol/kgSW
#   - DIC      :   Dissolved Inorganic Carbon in mol/kgSW
#   - NO3      :   Nitrate concentration in mol/kgSW
#   - SIOH4    :   Total Silicate in mol/kgSW
#
# Output : a data frame with the following colums
#    - CT      :   Conservative Temperature
#    - SA      :   Absolute Salinity
#
{
    # convert salinity
    SA <- sp2sa_chem (SP, TA=TA, DIC=DIC, NO3=NO3, SIOH4=SIOH4)
    # convert temperature
    CT <- gsw_CT_from_t (SA, T, P)
    
    return (as.data.frame(cbind(CT, SA)))
}


teos2eos_chem <- function (SA, CT, P=0, TA=2300e-6, DIC=2000e-6, NO3=0, SIOH4=0)
#
# This function converts Temperature and Salinity from the TEOS-10 to the EOS-80 standard.
# 
# Input parameters :
#   - SA       :   Absolute salinity in g/kg
#   - CT       :   Conservative Temperature in degrees C
#   - TA       :   Total Alkalinity in mol/kgSW
#   - DIC      :   Dissolved Inorganic Carbon in mol/kgSW
#   - NO3      :   Nitrate concentration in mol/kgSW
#   - SIOH4    :   Total Silicate in mol/kgSW
#
# Output : a data frame with the following colums
#    - T       :   in situ Temperature
#    - SP      :   Practical Salinity
#
{
    # convert temperature
    T <- gsw_t_from_CT (SA, CT, P)
    # convert salinity
    SP <- sa2sp_chem (SA, TA=TA, DIC=DIC, NO3=NO3, SIOH4=SIOH4)
    
    return (as.data.frame(cbind(T, SP)))
}


