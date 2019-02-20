# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 22:17:03 2019

@author: halle
"""

from calendarPlot import date_heatmap

import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
import datetime as dt
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon
#from astroplan import moon
plt.style.use(astropy_mpl_style)
import pandas as pd

### PARAMETERS ###
objectToObserve = SkyCoord.from_name('m33')
observerLocation = EarthLocation(lat=34.4242969*u.deg, lon=-119.8932443*u.deg, height=10*u.m)
utcoffset = -8*u.hour  # Pacific Daylight Time
OBJECT_MIN_ALT = 8 # degrees
MAX_MOON_ILLUMINATION = 0.1
MAX_MOON_ALT = 0 # degrees
MAX_SUN_ALT= -7 # degrees
startDate = dt.datetime(2019,1,1)
endDate = dt.datetime(2019,12,31)

def moon_illumination(thisTime):
    sun = get_sun(thisTime)
    moon = get_moon(thisTime)
    elongation = sun.separation(moon)
    i = np.arctan2(sun.distance*np.sin(elongation), moon.distance - sun.distance*np.cos(elongation))
    k = (1 + np.cos(i))/2.0
    return k

def angleBetweenPolarVectors(altAzi1,altAzi2):
    #delta phi = arccos(sin(phi_1)*sin(phi_2) + cos(phi_1)*cos(phi_2)*cos(lambda_2-lambda_1))
    deltaLambda = np.subtract(altAzi1.az.radian,altAzi2.az.radian)
    sines = np.multiply(np.sin(altAzi1.alt.radian),np.sin(altAzi2.alt.radian))
    cosines = np.multiply(np.cos(altAzi1.alt.radian),np.cos(altAzi2.alt.radian))
    return np.arccos(np.add(sines,np.multiply(cosines,np.cos(deltaLambda))))

def evaluateDay(thisDate):
    delta_midnight = np.linspace(-12, 12, 96)*u.hour
    tomorrow = thisDate + dt.timedelta(days=1)
    midnight = Time("{}-{}-{} 00:00:00".format(tomorrow.year,tomorrow.month,tomorrow.day)) - utcoffset
    moonIllum = moon_illumination(midnight) # assume constant illumination throughout the day
#    if (moonIllum > MAX_MOON_ILLUMINATION):
#        return 0
    allDay = midnight + delta_midnight
    allDayFrame = AltAz(obstime=allDay, location=observerLocation)
    objectAltAziAllDay = objectToObserve.transform_to(allDayFrame)
    objectAltAllDay = objectAltAziAllDay.alt
    sunAltAllDay = get_sun(allDay).transform_to(allDayFrame).alt
    moonAltAziAllDay = get_moon(allDay).transform_to(allDayFrame)
    moonAltAllDay = moonAltAziAllDay.alt
    #moonIllumAllDay = moon.moon_illumination(allDay)
    objectVisible = objectAltAllDay > (OBJECT_MIN_ALT * u.deg)
    nightTime = sunAltAllDay < (MAX_SUN_ALT * u.deg)
    moonLess = moonAltAllDay < (MAX_MOON_ALT * u.deg)
    goodHours = np.sum(objectVisible & nightTime & moonLess) * 0.25
    deltaAngles = angleBetweenPolarVectors(objectAltAziAllDay,moonAltAziAllDay) # Great circle angle between the moon and the object
#    fig1 = plt.figure(num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
#    ax1 = fig1.add_subplot(1,1,1)
#    ax1.plot(delta_midnight,objectAltAllDay,color='blue')
#    ax1.plot(delta_midnight,sunAltAllDay,color='red')
#    ax1.plot(delta_midnight,moonAltAllDay,color='black')
#    ax1.plot(delta_midnight,np.rad2deg(deltaAngles),color='green')
    return goodHours

def plotGoodness():
    dateData = []
    goodnessData = []
    dateData.append(startDate)
    numGoodHours = evaluateDay(startDate)
    goodnessData.append(numGoodHours)
    while dateData[-1] < endDate:
        dateData.append(dateData[-1] + dt.timedelta(days=1))
        print(dateData[-1].strftime("%B %d, %Y"))
        numGoodHours = evaluateDay(dateData[-1])
        goodnessData.append(numGoodHours)
    
    goodnessSeries = pd.Series(goodnessData)
    goodnessSeries.index =  pd.DatetimeIndex(dateData)
    # Create the figure. For the aspect ratio, one year is 7 days by 53 weeks.
    # We widen it further to account for the tick labels and color bar.
    figsize = plt.figaspect(7 / 60)
    fig = plt.figure(figsize=figsize)
    
    # Plot the heatmap with a color bar.
    ax = date_heatmap(goodnessSeries, edgecolor='black')
    #plt.colorbar(ticks=np.arange(0,1,0.2), pad=0.02)
    plt.colorbar()
    
    # Force the cells to be square. If this is set, the size of the color bar
    # may look weird compared to the size of the heatmap. That can be corrected
    # by the aspect ratio of the figure or scale of the color bar.
    ax.set_aspect('equal')
    annot = ax.annotate("", (1.25, 0.6),
                 xycoords="axes fraction", va="center", ha="center",
                 bbox=dict(boxstyle="round, pad=1", fc="w"))
    firstSunday = startDate - dt.timedelta(days=((startDate.weekday() + 1) % 7))
    
    def update_annot(thisX,thisY):
        thisTime = firstSunday + dt.timedelta(days=thisY + (thisX * 7))
        #text = "{:.2f}, {:.2f}".format(thisX,thisY)
        #thisMoonPhase = data.get(thisTime,0)
        numGoodHours = goodnessSeries.get(thisTime,0)
        text = "{}\n{:.2f} hours".format(thisTime.strftime("%B %d, %Y"),numGoodHours)
        annot.set_text(text)
        plt.draw()
    
    
    def hover(event):
    
        if not event.inaxes: 
                return
        
        xint = np.floor(event.xdata+0.5)
        yint = np.floor(event.ydata+0.5)
        update_annot(xint,yint)    
    
    fig.canvas.mpl_connect("motion_notify_event", hover)
    plt.show()
    
#start = time.time()
#print(evaluateDay(dt.datetime(2019,2,26)))
#end = time.time()
#print(end-start)
plotGoodness()