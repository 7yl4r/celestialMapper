from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np


class Measurement(object):
    '''
    This class represents a celestial object measurement,
    including its name and altitude above sea level.
    '''
    def __init__(self):
        self.altitude = None
        self.name = None
        self.recordMeasurement()
    
    def recordMeasurement(self):
        self.name = raw_input('enter name of object ')
        self.altitude = float(raw_input('enter observed altitude of '+self.name+' '))
        
class LOP(object):
    '''
    defines a Line Of Position (LOP) using lat/long and radius.
    '''
    def __init__(self,lat,long,rad):
        self.latitude = lat
        self.longitude = long
        self.radius = rad
        
    def draw(self, m, **kwargs):
        ''' 
        draws the LOP on the given map m.
        Given keyword arguments passed to matplotlib plot
        based on equi function found here: http://www.geophysique.be/2011/02/20/matplotlib-basemap-tutorial-09-drawing-circles/
        '''
        X = []
        Y = []
        for azimuth in range(0, 360):
            glon2, glat2, baz = self.__shoot(self.longitude, self.latitude, azimuth, self.radius)
            X.append(glon2)
            Y.append(glat2)
        X.append(X[0])
        Y.append(Y[0])
     
        #m.plot(X,Y,**kwargs) #Should work, but doesn't...
        X,Y = m(X,Y)
        plt.plot(X,Y,**kwargs)
        
        
    def __shoot(self, lon, lat, azimuth, maxdist=None):
        """Shooter Function
        Original javascript on http://williams.best.vwh.net/gccalc.htm
        Translated to python by Thomas Lecocq
        """
        glat1 = lat * np.pi / 180.
        glon1 = lon * np.pi / 180.
        s = maxdist / 1.852
        faz = azimuth * np.pi / 180.
     
        EPS= 0.00000000005
        if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
            alert("Only N-S courses are meaningful, starting at a pole!")
     
        a=6378.13/1.852
        f=1/298.257223563
        r = 1 - f
        tu = r * np.tan(glat1)
        sf = np.sin(faz)
        cf = np.cos(faz)
        if (cf==0):
            b=0.
        else:
            b=2. * np.arctan2 (tu, cf)
     
        cu = 1. / np.sqrt(1 + tu * tu)
        su = tu * cu
        sa = cu * sf
        c2a = 1 - sa * sa
        x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
        x = (x - 2.) / x
        c = 1. - x
        c = (x * x / 4. + 1.) / c
        d = (0.375 * x * x - 1.) * x
        tu = s / (r * a * c)
        y = tu
        c = y + 1
        while (np.abs (y - c) > EPS):
     
            sy = np.sin(y)
            cy = np.cos(y)
            cz = np.cos(b + y)
            e = 2. * cz * cz - 1.
            c = y
            x = e * cy
            y = e + e - 1.
            y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
                  d / 4. - cz) * sy * d + tu
     
        b = cu * cy * cf - su * sy
        c = r * np.sqrt(sa * sa + b * b)
        d = su * cy + cu * sy * cf
        glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
        c = cu * cy - su * sy * cf
        x = np.arctan2(sy * sf, c)
        c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
        d = ((e * cy * c + cz) * sy * c + y) * sa
        glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
     
        baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
     
        glon2 *= 180./np.pi
        glat2 *= 180./np.pi
        baz *= 180./np.pi
     
        return (glon2, glat2, baz)

def drawMap(*LOPs):
    ''' this function draws a map with given LOPs'''
    # create new figure, axes instances.
    fig=plt.figure()
    ax=fig.add_axes([0.1,0.1,0.8,0.8])
    # setup mercator map projection.
    map = Basemap(llcrnrlon=-100.,
        llcrnrlat=20.,
        urcrnrlon=20.,
        urcrnrlat=60.,\
        rsphere=(6378137.00,6356752.3142),\
        resolution='l',projection='merc',\
        lat_0=40.,lon_0=-20.,lat_ts=20.)
    map.drawcoastlines()
    # draw lat/lon grid lines every 30 degrees.
    map.drawmeridians(np.arange(0,360,30))
    map.drawparallels(np.arange(-90,90,30))
    for LOP in LOPs:
        LOP.draw(map)
    
def getLOP(measurement):
    ''' this function returns the LOP for given measurement'''
    #NOTE: currently this returns a dummy circle
    import random as rand
    centerlat = -30.231 + (rand.random()-.5)*20
    centerlon = 40.1412 + (rand.random()-.5)*20
    radius = 2000 + (rand.random()-.5)*1000
    return LOP(centerlon, centerlat, radius)
        
def main():
    ''' this function goes through the intercept_method'''
    print 'for 1st object'
    obj1 = Measurement()  
    print 'for 2nd object'
    obj2 = Measurement()
    
    LOP1 = getLOP(obj1)
    LOP2 = getLOP(obj2)
    drawMap(LOP1, LOP2)
    plt.show()
    print 'plotting complete.'
    
    
    
    

    
# this next bit makes main run if this file is called from cmd line:
if __name__ =='__main__':
    main()