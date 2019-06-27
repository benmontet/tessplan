import numpy as np
import matplotlib.pyplot as plt

# Query Exoplanet Orbit Database (exoplanets.org) for planet properties
from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astroplan import EclipsingSystem, Observer, FixedTarget, is_observable
from astroplan import (PrimaryEclipseConstraint, is_event_observable,
                       AtNightConstraint, AltitudeConstraint, LocalTimeConstraint)
import pytz

import urllib.parse as urlparse
import requests
from bs4 import BeautifulSoup
from tqdm import trange

from copy import deepcopy
from itertools import compress

__all__ = ['Planet', 'all_observable', 'make_local']

def make_local(times, timezone):
    tz = pytz.timezone('utc')
    tdt = times.to_datetime()
    outarr = np.array([])

    for i in range(len(times)):
        inplace = tz.localize(tdt[i])
        val = inplace.astimezone(pytz.timezone(timezone))
        outarr = np.append(outarr, val.strftime('%Y-%m-%d %H:%M:%S'))
    return outarr
 

def all_observable(date, run_length, **kwargs):
    
    all_obs = np.array([])
    
    url = 'https://exofop.ipac.caltech.edu/tess/view_toi.php'
    page = requests.get(url)
    
    soup = BeautifulSoup(page.content, 'html.parser')
    table = list(soup.children)[9]
    
    toitable = list(table.children)[0].split('{')[108:-49]
    
    for j in trange(len(toitable)):
        listone = toitable[j].split('}')[0].split(',')
        
            
        truefalse = [':' in listone[i] for i in range(len(listone))]
        listtwo = deepcopy(listone)
        for i in range(len(listone)-1, -1, -1):
            if truefalse[i] == False:
                listtwo[i-1] = listtwo[i-1] + ',' + listtwo[i]
        list_all = list(compress(listtwo, truefalse))
        
        
        keys = [list_all[i].split(':')[0].lstrip() for i in range(len(list_all))]
        values = [list_all[i].split(':')[1].lstrip() for i in range(len(list_all))]
        dd ={ keys[i] : values[i] for i in range(0, len(keys) ) }
        
        try:


            goodp = Planet(t0=float(dd['epoch'])-2457000, period=float(dd['period']), perr=float(dd['period_e']), duration=float(dd['duration']),
                           ra=float(dd['ra']), dec=float(dd['dec']), startdate=date,run_length=run_length, **kwargs)

            goodp.depth = float(dd['depth_ppm'])
            goodp.disposition = dd['disposition']
            goodp.disposition_tfop = dd['tfopwg_disposition']
            try:
                goodp.pl_rad = float(dd['planet_radius'])
            except:
                goodp.pl_rad = None
            try:
                goodp.st_rad = float(dd['stellar_radius'])
            except:
                goodp.st_rad = None
            try:
                goodp.logg = float(dd['stellar_logg'])
            except:
                goodp.logg = None
            goodp.teff = float(dd['stellar_teff'])
            goodp.tmag = float(dd['tess_mag'])
            goodp.tic = int(dd['ticid'][1:-1])
            goodp.toi = dd['toi']
            goodp.sg1a = int(dd['sg1a_priority'])
            goodp.sg1b = int(dd['sg1b_priority'])
            goodp.sg2 = int(dd['sg2_priority'])
            goodp.sg3 = int(dd['sg3_priority'])
            goodp.sg4 = int(dd['sg4_priority'])
            goodp.sg5 = int(dd['sg5_priority'])


            if goodp.obs_mid_times_utc is not None:
                if len(goodp.obs_mid_times_utc) > 0:
                    all_obs = np.append(all_obs, goodp)
        except:
            pass

    return all_obs


class Planet(object):
    
    def __init__(self, t0=None, period=None, perr=None, duration=None, loc = 'Siding Spring Observatory', timezone='Australia/NSW',
                ra=None, dec=None, startdate=None, starttime='0:00', run_length=180, el_limit=30, toi=None):
    
        
        if toi is not None:
            self.set_params_from_toi(toi)

                    
        else:
            
            self.epoch = Time(t0+2457000, format='jd')
            self.period = period * u.d
            self.period_err = perr * u.d

            self.duration = duration / 24 * u.d
            self.coords = [ra*u.deg, dec*u.deg]
        
        self.alt_limit = el_limit
        
        self.observatory = loc
        self.timezone = timezone
        
        
        
        self.obs_mid_times_utc = None
        self.obs_mid_uncerts = None
        self.obs_mid_times_local = None
        
        
        
        coord = SkyCoord(ra=self.coords[0], dec=self.coords[1])
        target = FixedTarget(coord, name='Target')


        SSO = Observer.at_site(self.observatory, timezone=self.timezone)
        


        planet = EclipsingSystem(primary_eclipse_time=self.epoch, orbital_period=self.period,
                                   duration=self.duration)
        
        starttime = startdate + ' ' + starttime
        self.start = Time(starttime)
        self.run_length = run_length * u.d

        n_transits = np.ceil(self.run_length / self.period)  
        

        self.all_transit_times = planet.next_primary_eclipse_time(self.start, n_eclipses=n_transits)
        
        diff = self.all_transit_times - self.start - self.run_length
        real_trans = diff.value < 0
        if np.sum(real_trans) < 0.5:
            self.all_transit_times = None
            return
        else:
            self.all_transit_times = self.all_transit_times[real_trans]
            

        self.uncert_vals = (self.all_transit_times - self.epoch)*self.period_err/self.period * 1440

        constraints = [AtNightConstraint.twilight_nautical(),
                AltitudeConstraint(min=self.alt_limit*u.deg)]

        obs_mid = is_event_observable(constraints, SSO, target, times=self.all_transit_times)[0]
        obs_start = is_event_observable(constraints, SSO, target, times=self.all_transit_times-0.5*self.duration)[0]
        obs_end = is_event_observable(constraints, SSO, target, times=self.all_transit_times+0.5*self.duration)[0]
        

        if np.sum(obs_mid*obs_start*obs_end) > 0:
            
            self.obs_airmass = np.zeros((np.sum(obs_mid*obs_start*obs_end), 3))

            self.obs_mid_times_utc = self.all_transit_times[obs_mid*obs_start*obs_end]
            self.obs_mid_uncerts   = self.uncert_vals[obs_mid*obs_start*obs_end]
            
            for i in range(np.sum(obs_mid*obs_start*obs_end)):
                obs_airmass_mid = SSO.altaz(self.obs_mid_times_utc[i], coord).secz 
                obs_airmass_start = SSO.altaz(self.obs_mid_times_utc[i]-self.duration/2, coord).secz 
                obs_airmass_end = SSO.altaz(self.obs_mid_times_utc[i]+self.duration/2, coord).secz 

                self.obs_airmass[i] = [obs_airmass_start, obs_airmass_mid, obs_airmass_end]


            # NEED A BARYCENTRIC CORRECTION

            tz = pytz.timezone('utc')


            dtime = self.obs_mid_times_utc.to_datetime()

            self.obs_mid_times_local = np.array([])

            for i in range(len(self.obs_mid_times_utc)):
                inoz = tz.localize(dtime[i])
                indiv_uncert = self.obs_mid_uncerts[i]
                val = inoz.astimezone(pytz.timezone(self.timezone))
                self.obs_mid_times_local = np.append(self.obs_mid_times_local, val.strftime('%Y-%m-%d %H:%M:%S'))


    
    def set_params_from_toi(self, true_toi):
        true_toi = str(true_toi)
        url = 'https://exofop.ipac.caltech.edu/tess/view_toi.php'
        page = requests.get(url)

        soup = BeautifulSoup(page.content, 'html.parser')
        table = list(soup.children)[9]

        toitable = list(table.children)[0].split('{')[108:-49]
        

        for j in range(len(toitable)):
            toi = toitable[j].split('}')[0].split(',')[1].split(':')[1].lstrip()
            
            if toi == true_toi:

                listone = toitable[j].split('}')[0].split(',')
                truefalse = [':' in listone[i] for i in range(len(listone))]
                listtwo = deepcopy(listone)
                for i in range(len(listone)-1, -1, -1):
                    if truefalse[i] == False:
                        listtwo[i-1] = listtwo[i-1] + ',' + listtwo[i]
                list_all = list(compress(listtwo, truefalse))


                keys = [list_all[i].split(':')[0].lstrip() for i in range(len(list_all))]
                values = [list_all[i].split(':')[1].lstrip() for i in range(len(list_all))]
                dd ={ keys[i] : values[i] for i in range(0, len(keys) ) }


                t0 = float(dd['epoch'])
                self.epoch = Time(t0, format='jd')
                try:
                    self.period = float(dd['period']) * u.d
                    self.period_err = float(dd['period_e']) * u.d
                except:
                    self.period = None
                    self.period_err = None
                self.duration = float(dd['duration']) * u.d
                ra = float(dd['ra'])
                dec = float(dd['dec'])
                self.coords = [ra*u.deg, dec*u.deg]
                self.depth = float(dd['depth_ppm'])
                self.disposition = dd['disposition']
                self.disposition_tfop = dd['tfopwg_disposition']
                try:
                    self.pl_rad = float(dd['planet_radius'])
                except:
                    self.pl_rad = None
                try:
                    self.st_rad = float(dd['stellar_radius'])
                except:
                    self.st_rad = None
                try:
                    self.logg = float(dd['stellar_logg'])
                except:
                    self.logg = None
                self.teff = float(dd['stellar_teff'])
                self.tmag = float(dd['tess_mag'])
                self.tic = int(dd['ticid'][1:-1])
                self.toi = dd['toi']
                self.sg1a = int(dd['sg1a_priority'])
                self.sg1b = int(dd['sg1b_priority'])
                self.sg2 = int(dd['sg2_priority'])
                self.sg3 = int(dd['sg3_priority'])
                self.sg4 = int(dd['sg4_priority'])
                self.sg5 = int(dd['sg5_priority'])
            
                return
        print('TOI not found!')
