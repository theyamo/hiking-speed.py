#!/usr/bin/python3
import math
from lxml import etree
import datetime
import optparse

def tobler_function(slope):
    return 6 * math.exp(-3.5 * abs(math.tan(math.radians(slope)) + 0.05))

def calc_slope(elev, distance):
    return math.atan(elev / distance) * 360 / (2 * math.pi)

def calc_distance(lat1, lon1, lat2, lon2):
    """
    Calculates geodetic distance between two points specified by latitude/longitude using
    Vincenty inverse formula for ellipsoids
    
    @param lat1: latitude first point, in decimal degrees
    @param lon1: longitude first point, in decimal degrees
    @param lat2: latitude second point, in decimal degrees
    @param lon2: longitude second point, in decimal degrees
    @return: distance in metres between points
    """

    WGS84_a = 6378137
    WGS84_b = 6356752.314245
    WGS84_f = 1/298.257223563  # WGS-84 ellipsoid params
    L = math.radians(lon2 - lon1)
    U1 = math.atan((1-WGS84_f) * math.tan(math.radians(lat1)))
    U2 = math.atan((1-WGS84_f) * math.tan(math.radians(lat2)))
    sinU1 = math.sin(U1)
    cosU1 = math.cos(U1)
    sinU2 = math.sin(U2)
    cosU2 = math.cos(U2)
    lambda_ = L
    iterlimit = 20
    
    while True:
        sinLambda = math.sin(lambda_)
        cosLambda = math.cos(lambda_)
        sinSigma = math.sqrt((cosU2*sinLambda) * (cosU2*sinLambda) +
            (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda))
        
        if sinSigma==0.0:
            return 0.0 # co-incident points
        
        cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda
        sigma = math.atan2(sinSigma, cosSigma)
        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
        cosSqAlpha = 1 - sinAlpha*sinAlpha
        cos2SigmaM = cosSigma - 2*sinU1*sinU2/cosSqAlpha

        # equatorial line: cosSqAlpha=0 (6)
        if math.isnan(cos2SigmaM):
            cos2SigmaM = 0.0
        
        C = WGS84_f/16 * cosSqAlpha * (4 + WGS84_f * (4 - 3 * cosSqAlpha))
        lambdaP = lambda_
        lambda_ = L + (1-C) * WGS84_f * sinAlpha * (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1 + 2 * cos2SigmaM * cos2SigmaM)))

        iterlimit -= 1
        if not abs(lambda_-lambdaP) > 1e-12 and iterlimit > 0:
            break

    if iterlimit==0:
        return None  # NB: returned NAN // formula failed to converge
    
    uSq = cosSqAlpha * (WGS84_a*WGS84_a - WGS84_b*WGS84_b) / (WGS84_b*WGS84_b)
    A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
    B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
    deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-
        B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)))
    s = WGS84_b*A*(sigma-deltaSigma)
    
    s = round( s * 1000 ) / 1000 # round to 1mm precision
    return s

def timestamp(time):
    seconds = round(time) % 60
    mins = math.trunc(time / 60) % 60
    hours = math.trunc(time / 60 / 60) % 24
    days = math.trunc(time / 60 / 60 / 24)
    return '%02d:%02d:%02d:%02d' % (days, hours, mins, seconds)

def format_utc_timestamp(time):
    seconds = math.trunc(time % 60)
    mins = math.trunc(time / 60) % 60
    hours = math.trunc(time / 60 / 60) % 24
    days = math.trunc(time / 60 / 60 / 24)
    dt = datetime.datetime(1970, 1, days + 1, hours, mins, seconds)
    return dt.strftime('%Y-%m-%dT%H:%M:%S.000Z')

def go(infn, outfn):
    ns = {'gpx': 'http://www.topografix.com/GPX/1/1',
          'gpxx': '"http://www.garmin.com/xmlschemas/GpxExtensions/v3' }

    root = etree.parse(infn).getroot()

    prevlat = prevlon = prevele = 0.0
    total_distance = 0
    first_iter = True
    total_time = 0
    point = 0

    for trkpt in root.findall('.//*/gpx:trkpt', ns):
        lat = float(trkpt.get('lat'))
        lon = float(trkpt.get('lon'))
        dist = 0

        if not first_iter:
            dist = calc_distance(prevlat, prevlon, lat, lon)

        first_iter = False
            
        if not math.isnan(dist):
            total_distance += dist

        ele = float('nan')
        for elept in trkpt.findall('./gpx:ele', ns):
            ele = float(elept.text)

        slope = 0.0
        if not math.isnan(ele) and dist > 0:
            slope = calc_slope(ele - prevele, dist)

        velocity = 0.0
        if not math.isnan(slope):
            velocity = round(tobler_function(slope), 2)

        # meters per second
        mps = velocity * 1000 / 60 / 60
        total_time += dist / mps

        if not outfn:
            print('#%d lat=%.4f lon=%.4f dist=%.3f tot=%.3f ele=%.2f slope=%.3f vel=%.2f t=%s'
                  % (point,
                     round(lat, 4),
                     round(lon, 4),
                     dist,
                     round(total_distance, 3),
                     ele,
                     round(slope, 3),
                     round(velocity, 3),
                     timestamp(total_time)))

        # get time or override if exists
        dt = format_utc_timestamp(round(total_time))

        node = trkpt.findall('time', ns)
        if len(node) > 1:
            raise Exception('cannot parse xml')
        elif len(node) == 1:
            node = node[0]
        elif len(node) == 0:
            node = etree.SubElement(trkpt, 'time')
        node.text = dt

        node = trkpt.findall('speed', ns)
        if len(node) > 1:
            raise Exception('cannot parse xml')
        elif len(node) == 1:
            node = node[0]
        elif len(node) == 0:
            node = etree.SubElement(trkpt, 'speed')
        node.text = str(velocity)

        prevlat = lat
        prevlon = lon
        prevele = ele
        point += 1

    if not outfn:
        print('total distance = ', round(total_distance, 1))

    # write xml
    if(outfn):
        tree = etree.ElementTree(root)
        tree.write(open(outfn, 'wb'))

if __name__ == '__main__':
    parser = optparse.OptionParser('usage: %prog [options]')
    parser.add_option('-o', '--output', dest = 'outfn', default = None, type = 'string', help = 'output to gpx file')
    (options, args) = parser.parse_args()

    infn = ''
    outfn = None

    if(len(args) < 1):
        print('error: filename required')
        exit(-1)
    
    infn = args[0]
    outfn = options.outfn

    go(infn, outfn)


