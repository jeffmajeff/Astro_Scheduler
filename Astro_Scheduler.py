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
from matplotlib import gridspec
import matplotlib.font_manager as fm

### PARAMETERS ###
objectToObserve = SkyCoord.from_name('m33')
observerLocation = EarthLocation(lat=34.4242969*u.deg, lon=-119.8932443*u.deg, height=10*u.m)
utcoffset = -8*u.hour  # Pacific Daylight Time
OBJECT_MIN_ALT = 8 # degrees
MAX_MOON_ILLUMINATION = 0.1
MAX_MOON_ALT = 0 # degrees
MAX_SUN_ALT= -7 # degrees
startDate = dt.datetime(2019,1,1)
endDate = dt.datetime(2019,1,31)

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
    #deltaAngles = np.rad2deg(angleBetweenPolarVectors(objectAltAziAllDay,moonAltAziAllDay)) # Great circle angle between the moon and the object
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
    #figsize = plt.figaspect(7 / 60)
    fig = plt.figure(figsize=(15,3))
    gs = gridspec.GridSpec(1, 16)
    
    # Plot the heatmap with a color bar.
    ax1 = plt.subplot(gs[0:12])
    ax2 = plt.subplot(gs[12])
    ax2.axis('off')
    ax3 = plt.subplot(gs[13:])
    ax1 = date_heatmap(goodnessSeries, edgecolor='black',ax=ax1)
    #plt.colorbar(ticks=np.arange(0,1,0.2), pad=0.02)
    plt.colorbar()
    
    # Force the cells to be square. If this is set, the size of the color bar
    # may look weird compared to the size of the heatmap. That can be corrected
    # by the aspect ratio of the figure or scale of the color bar.
    ax1.set_aspect('equal')
    annot = ax2.annotate("", (0.1, 0.5),
                 xycoords="axes fraction", va="center", ha="center",
                 bbox=dict(boxstyle="round, pad=1", fc="w"), fontsize=8)
    firstSunday = startDate - dt.timedelta(days=((startDate.weekday() + 1) % 7))
    
    font = fm.FontProperties(size=8)
    for label in (ax3.get_xticklabels() + ax3.get_yticklabels()):
        label.set_fontproperties(font)
    
    def update_annot(thisX,thisY):
        thisTime = firstSunday + dt.timedelta(days=thisY + (thisX * 7))
        #text = "{:.2f}, {:.2f}".format(thisX,thisY)
        #thisMoonPhase = data.get(thisTime,0)
        numGoodHours = goodnessSeries.get(thisTime,0)
        text = "{}\n{:.2f} hours".format(thisTime.strftime("%B %d, %Y"),numGoodHours)
        annot.set_text(text)
        plt.draw()
    
    
    def hover(event):
    
        if event.inaxes != ax1: 
                return
        
        xint = np.floor(event.xdata+0.5)
        yint = np.floor(event.ydata+0.5)
        update_annot(xint,yint)
        
    def onclick(event):
        #if not event.inaxes: 
        if event.inaxes != ax1:
                return
        xint = np.floor(event.xdata+0.5)
        yint = np.floor(event.ydata+0.5)
        delta_midnight = np.linspace(-12, 12, 96)*u.hour
        tomorrow = firstSunday + dt.timedelta(days=yint + (xint * 7)) + dt.timedelta(days=1)
        midnight = Time("{}-{}-{} 00:00:00".format(tomorrow.year,tomorrow.month,tomorrow.day)) - utcoffset

        allDay = midnight + delta_midnight
        allDayFrame = AltAz(obstime=allDay, location=observerLocation)
        objectAltAziAllDay = objectToObserve.transform_to(allDayFrame)
        objectAltAllDay = objectAltAziAllDay.alt
        sunAltAllDay = get_sun(allDay).transform_to(allDayFrame).alt
        moonAltAziAllDay = get_moon(allDay).transform_to(allDayFrame)
        moonAltAllDay = moonAltAziAllDay.alt
        deltaAngles = np.rad2deg(angleBetweenPolarVectors(objectAltAziAllDay,moonAltAziAllDay)) # Great circle angle between the moon and the object
        ax3.clear()
        ax3.plot(delta_midnight,objectAltAllDay,color='blue')
        ax3.plot(delta_midnight,sunAltAllDay,color='red')
        ax3.plot(delta_midnight,moonAltAllDay,color='black')
        ax3.plot(delta_midnight,deltaAngles,color='green')
        ax3.set_title((firstSunday + dt.timedelta(days=yint + (xint * 7))).strftime("%B %d, %Y"), fontsize=8)
        ax3.set_ylabel('degrees', fontsize=8)
        ax3.set_xlabel('hours around midnight', fontsize=8)
        ax3.set_ylim(0,np.max([90,np.max(deltaAngles)]))
        
        

    fig.canvas.mpl_connect('button_press_event', onclick)
    fig.canvas.mpl_connect("motion_notify_event", hover)
    plt.show()
    
#start = time.time()
#print(evaluateDay(dt.datetime(2019,2,26)))
#end = time.time()
#print(end-start)
plotGoodness()
