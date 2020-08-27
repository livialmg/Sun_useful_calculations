#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 10:54:51 2020

This scripts will implement some equations for solar position calculations

@author: Livia Gava, PGMET-INPE
"""
import numpy as np

def calc_declination(dayJ, convert = None):
    """
    INPUT:
    dayJ: Julian day
    OUTPUT:
    delta: Solar declination in radians (to convert to degrees, convert=True)
    """
    #      Parametros geometricos
    pi = 3.141592653589793
    convgr = 180. / pi
    #Calculo da declinacao (Paltridge e Platt 1976, baseado em Spencer 1971).
    #this equation result is in radians
    teta = (2. * pi * (dayJ - 1.) / 365.)       
    delta = (.006918 - (.399912 * np.cos(teta))
    + (.070257 * np.sin(teta)) - (.006758 * np.cos(2. * teta))
    + (.000907 * np.sin(2. * teta)) - (.002697 * np.cos(3. * teta))
    + (.00148 * np.sin(3. * teta)))
    if convert:
        delta = delta * convgr
    return delta
    
def calc_eqtime(dayJ, convert = None):
    """
    INPUT:
    dayJ: Julian day
    OUTPUT:
    eqtime: time correction from time eqaution in radians
    (to convert to minutes, convert = 'MIN', to hour convert = 'HOUR')
    """
    #      Parametros geometricos
    pi = 3.141592653589793

    teta = (2. * pi * (dayJ - 1.) / 365.)   
    #Equacao do tempo (Spencer 1971 in Paltridge e Platt 1976).
    #this equation result is in radians (LIOU)
    eqtime = (0.000075 + 0.001868 * np.cos(teta) 
    - 0.032077 * np.sin(teta) - 0.014615 * np.cos(2. * teta)
    - 0.040849 * np.sin(2. * teta))
    if convert == 'MIN':
        eqtime = eqtime * 229.18 #eqtime in minutes
    elif convert == 'HOUR':
        eqtime = eqtime * (180. / (pi * 15.)) #eqtime in hours
    return eqtime

def true_solar_time(long, hour, dayJ, localtime = None, timezone = None, result = None):
    """ 
    INPUT:
    long: longitude in degrees (positive to the east of the Prime Meridian)
    hour: hours in minutes
    if hour in local time, set localtime = True and gives the correspondent time zone (in hours from UTC)
    if hour in UTC localtime = None 
    dayJ: Julian day
    OUTPUT:
    tst: True Solar Time (to results in hour, set result = 'hour', else result will be in minutes)
    """
    if localtime:
        if result == 'hour':
            eqtime = calc_eqtime(dayJ, convert = 'HOUR')
            tst = (hour / 60.) + eqtime + ((4. * long) / 60.) - timezone
        else:
            eqtime = calc_eqtime(dayJ, convert = 'MIN')
            tst = hour + eqtime + (4. * long) - (60. * timezone)
    else:
        if result == 'hour':
            eqtime = calc_eqtime(dayJ, convert = 'HOUR')
            tst = (hour / 60.) + eqtime + (long / 15.)
        else:
            eqtime = calc_eqtime(dayJ, convert = 'MIN')
            tst = hour + eqtime + ((long / 15.) * 60.)
    return tst

def calc_hour_angle(true_solar_time, in_minute = None, in_degrees = None):
    """
    INPUT:
    true_solar_time: the solar time in hour, if tst in minutes set in_minute = True
    OUTPUT:
    ha: hour angle in radians (morning is negative), to hour angle in degree, set in_degrees = True
    """
    pi = 3.141592653589793

    if in_minute:
        ha = (true_solar_time / 4.) - 180. 
        if in_degrees:
            ha = ha
        else:
            ha = ha * (pi / 180.)
    else:
        ha = (true_solar_time - 12.) * 15.
        if in_degrees:
            ha = ha
        else:
            ha = ha * (pi / 180.)
    return ha

def calc_cosZ(lat, lon, dayJ, hour, hour_in_minutes = None):
    """
    INPUT:
    lat: latitude in degrees
    lon: longitude in degrees
    dayJ: Julian day
    hour: hour in UTC (if passed in minutes (1 to 1440) set hour_in_minutes = True)
    OUTPUT:
    mu0: solar zenith angle cosine
    """
    # Parametros geometricos
    pi = 3.141592653589793
    convgr = pi / 180. 
    
    ylat = lat * convgr #converts degrees to radians
    
    if hour_in_minutes:
        delta = calc_declination(dayJ)
        tst = true_solar_time(lon, hour, dayJ) #true solar time in minutes
        hour_angle = calc_hour_angle(tst, in_minute = True)
        #Calcula Angulo Zenital:
        #(Varejao-Silva e Ceballos, 1980)
        mu0 = np.cos(delta) * np.cos(ylat) * np.cos(hour_angle) + np.sin(delta) * np.sin(ylat)
    else:
        hour_in_minutes = hour * 60
        delta = calc_declination(dayJ)
        tst = true_solar_time(lon, hour_in_minutes, dayJ)
        hour_angle = calc_hour_angle(tst, in_minute = True)
        #Calcula Angulo Zenital:
        #(Varejao-Silva e Ceballos, 1980)
        mu0 = np.cos(delta) * np.cos(ylat) * np.cos(hour_angle) + np.sin(delta) * np.sin(ylat)
    return mu0

def calc_azimuth(lat, lon, dayJ, hour, hour_in_minutes = None):
    """
    INPUT:
    lat: latitude in degrees
    lon: longitude in degrees
    dayJ: Julian day
    hour: hour in UTC (if passed in minutes (1 to 1440) set hour_in_minutes = True)
    OUTPUT:
    azim:  solar azimuth (degrees clockwise from north)
    """
    # Parametros geometricos
    pi = 3.141592653589793
    convgr = pi / 180. 
    
    ylat = lat * convgr #converts degrees to radians
    
    if hour_in_minutes:
        delta = calc_declination(dayJ) #radians
        mu0 = calc_cosZ(lat, lon, dayJ, hour, hour_in_minutes = True)
        Z = np.arccos(mu0)
        cos_azim = (((np.sin(ylat) * mu0) - np.sin(delta)) / (np.sin(Z) * np.cos(ylat)))
        azim = np.arccos(cos_azim)
    else:
        delta = calc_declination(dayJ) #radians
        mu0 = calc_cosZ(lat, lon, dayJ, hour)
        Z = np.arccos(mu0)
        cos_azim = (((np.sin(ylat) * mu0) - np.sin(delta)) / (np.sin(Z) * np.cos(ylat)))
        azim = np.arccos(cos_azim)
    return azim
        
def distance_sun_earth(dayJ):
    """
    INPUT:
    dayJ: Julian day
    OUTPUT:
    dist: distance Sun-Earth for a given day
    """
    med_dist = 149600000 #mean sun-earth distance in kilometers
    pi = 3.141592653589793
    #from Spencer.
    teta = (2. * pi * (dayJ - 1.) / 365.) 
    dist = med_dist / np.sqrt(1.000110 + 0.034221 * np.cos(teta) + 0.001280 * np.sin(teta) +
                              0.000719 * np.cos(2. * teta) + 0.000077 * np.sin(2 * teta))
    return dist

def daylengh(lat, dayJ):
    """
    INPUT:
    lat: latitude in degrees
    dayJ: Julian day
    OUTPUT:
    N: Max number of minutes of day light
    """
    pi = 3.141592653589793
    convgr = pi / 180.
    
    ylat = lat * convgr #converts latitude from degrees to radians
    zenith_corrected = 90.833 * (pi / 180.)
    delta = calc_declination(dayJ)    
    #Calculates the day lenght in minutes
    N = (24. / pi) * np.arccos((np.cos(zenith_corrected) / 
                    (np.cos(ylat) * np.cos(delta))) - (np.tan(ylat) * np.tan(delta))) * 60.    
    return N


def calc_sunrise(lat, lon, dayJ, result_in_hour = None):
    """
    INPUT:
    lat: latitude in degrees
    long: longitude in degrees
    dayJ: Julian day
    OUTPUT:
    sunrise: UTC time of sunrise in minutes (to output in HH:MM format, result_in_hour = True)
    ha_degree: hour angle of sunrise in degree
    """
    pi = 3.141592653589793
    ylat = lat * (pi / 180.)
    delta = calc_declination(dayJ)
    eqtime = calc_eqtime(dayJ, convert = 'MIN')
    zenith_corrected = 90.833 * (pi / 180.)
    ha = np.arccos((np.cos(zenith_corrected) / 
                    (np.cos(ylat) * np.cos(delta))) - (np.tan(ylat) * np.tan(delta))) #in radians
    ha_degree = ha * (180. / pi)
    sunrise = 720 - (4 * (lon + ha_degree)) - eqtime
    if result_in_hour:
        result = convert_minutes_to_hour(sunrise, hour_format = True)
        return result, ha_degree
    else:
        return sunrise, ha_degree

def calc_sunset(lat, lon, dayJ, result_in_hour = None):
    """
    INPUT:
    lat: latitude in degrees
    long: longitude in degrees
    dayJ: Julian day
    OUTPUT:
    sunset: UTC time of sunset in minutes (to output in HH:MM format, result_in_hour = True)
    ha_degree: hour angle of sunset in degrees
    """
    pi = 3.141592653589793
    ylat = lat * (pi / 180.)
    delta = calc_declination(dayJ)
    eqtime = calc_eqtime(dayJ, convert = 'MIN')
    zenith_corrected = 90.833 * (pi / 180.)
    ha = - (np.arccos((np.cos(zenith_corrected) / 
                    (np.cos(ylat) * np.cos(delta))) - (np.tan(ylat) * np.tan(delta)))) #in radians
    ha_degree = ha * (180. / pi)
    sunset = 720 - (4 * (lon + ha_degree)) - eqtime
    if result_in_hour:
        result = convert_minutes_to_hour(sunset, hour_format = True)
        return result, ha_degree
    else:
        return sunset, ha_degree

def solar_noon(lon, dayJ):
    """
    long: longitude in degrees
    dayJ: Julian day
    OUTPUT:
    solar_noon: Solar noon for a given location
    """
    eqtime = calc_eqtime(dayJ, convert = 'MIN')
    solar_noon = 720 - 4*lon - eqtime
    return solar_noon
    
def convert_minutes_to_hour(hour_in_minutes, hour_format = None):
    """
    INPUT:
    hour_in_minutes: hour given in minutes (1 to 1440)
    OUTPUT:
    for hour in HH:MM format, set hour_format = True, otherwise returns a tuple with the hour, 
    and minutes as (hour, minutes).
    """
    hour = int(hour_in_minutes // 60)
    minute = round(hour_in_minutes % 60)
    
    if hour_format:
        return "%02d:%02d" % (hour, minute)
    else:
        return hour, minute