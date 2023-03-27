"""
Module that contains the ARCSILISS LISS III Sensor class.
"""
############################################################################
#  arcsisensorLISS III.py
#
#  Copyright 2013 ARCSI.
#
#  ARCSI: 'Atmospheric and Radiometric Correction of Satellite Imagery'
#
#  ARCSI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ARCSI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with ARCSI.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Purpose:  A class for read the LISS III sensor header file and applying
#           the pre-processing operations within ARCSI to the LISS III 8
#           datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 05/07/2013
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

# Import updated print function into python 2.7
from __future__ import print_function
# Import updated division operator into python 2.7
from __future__ import division
# import abstract base class stuff
from .arcsisensor import ARCSIAbstractSensor
# Import the ARCSI exception class
from .arcsiexception import ARCSIException
# Import the ARCSI utilities class
from .arcsiutils import ARCSIUtils
# Import the datetime module
import datetime
# Import the GDAL/OGR spatial reference library
from osgeo import osr
from osgeo import ogr
# Import OS path module for manipulating the file system
import os.path
# Import the RSGISLib Module.
import rsgislib
# Import the RSGISLib Image Calibration Module.
import rsgislib.imagecalibration
# Import the RSGISLib Image Utilities Module.
import rsgislib.imageutils
# Import the RSGISLib Image Calculations Module.
import rsgislib.imagecalc
# Import the RSGISLib segmentation Module
import rsgislib.segmentation
# Import the RSGISLib segmentation Module
import rsgislib.segmentation.segutils
# Import the RSGISLib Raster GIS Module
import rsgislib.rastergis
# Import the collections module
import collections
# Import the py6s module for running 6S from python.
import Py6S
# Import the python maths library
import math
# Import the RIOS RAT library
from rios import rat
# Import the GDAL python library
import osgeo.gdal as gdal
# Import the scipy optimisation library - used for finding AOD values form the imagery.
from scipy.optimize import minimize
# Import the numpy library
import numpy
# Import JSON module
import json
# Import the shutil module
import shutil
# Import the solar angle tools from RSGISLib
import rsgislib.imagecalibration.solarangles
import rios.fileinfo
#import sunpy module for earth sun distance calculation
from sunpy import coordinates



class ARCSILISS_Sensor(ARCSIAbstractSensor):
    """
    A class which represents the LISS III 8 sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "L3"
        
        self.band3File = ""
        self.band4File = ""
        self.band5File = ""
        self.band6File = ""
       
        self.row = 0
        self.path = 0
        
        self.b3MinRad = 0.0
        self.b3MaxRad = 0.0 
        self.b4MinRad = 0.0
        self.b4MaxRad = 0.0
        self.b5MinRad = 0.0
        self.b5MaxRad = 0.0
        self.b6MinRad = 0.0
        self.b6MaxRad = 0.0
        
        
        
        
        self.b3RadMulti = 0.0
        self.b3CalMax = 0
        self.b4RadMulti = 0.0
        self.b4CalMax = 0
        self.b5RadMulti = 0.0
        self.b5CalMax = 0
        self.b6RadMulti = 0.0
        self.b6CalMax = 0
        

       
        self.b3RadAdd = 0.0
        self.b3MaxRad = 0.0
        self.b4RadAdd = 0.0
        self.b4MaxRad = 0.0
        self.b5RadAdd = 0.0
        self.b5MaxRad = 0.0
        self.b6RadAdd = 0.0
        self.b6MaxRad = 0.0
      
                    
        self.b3CalMin = 0.0
        self.b3CalMax = 0.0
        self.b4CalMin = 0.0
        self.b4CalMax = 0.0
        self.b5CalMin = 0.0
        self.b5CalMax = 0.0
        self.b6CalMin = 0.0
        self.b6CalMax = 0.0
             
        self.sensorID = ""
        self.spacecraftID = ""        
        self.earthSunDistance = 0.0       
        self.gridCellSizeRefl = 0.0
        self.satAltitude = 0.0
       
    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the LISS III MTL header files
        """
        try:
            if not self.userSpInputImage is None:
                raise ARCSIException("LISS III sensor cannot accept a user specified image file - only the images in the header file will be used.")
            self.headerFileName = os.path.split(inputHeader)[1]
            arcsiUtils = ARCSIUtils()            
            print("Reading header file")
            hFile = open(inputHeader, 'r')
            headerParams = dict()
            for line in hFile:
                line = line.strip()
                if line:
                    lineVals = line.split('=')
                    if len(lineVals) == 2:
                        if (lineVals[0].strip() != "JpegNoRows"):
                            headerParams[lineVals[0].strip()] = lineVals[1].strip().replace('"','')
            hFile.close()
            print("Extracting Header Values")
            # Get the sensor info.
            if ((headerParams["SatID"] == "IRS-R2") and (headerParams["Sensor"] == "L3")):
                self.sensor = "L3"
            else:
                raise ARCSIException("Do no recognise the spacecraft and sensor or combination.")
            self.sensorID = headerParams["Sensor"]
            self.spacecraftID = headerParams["SatID"]
            # Get row/path
            self.row = int(headerParams["Row"])
            self.path = int(headerParams["Path"])
            # Get date and time of the acquisition
            #acData1 = datetime.datetime.strptime(headerParams["DateOfPass"], '%d-%b-%Y')
            acData = datetime.date.isoformat(datetime.datetime.strptime(headerParams["DateOfPass"], '%d-%b-%Y')).split('-')
            acTime = headerParams["SceneCenterTime"][12:-1].split(':')
            secsTime = acTime[2].split('.')
            self.acquisitionTime = datetime.datetime(int(acData[0]), int(acData[1]), int(acData[2]), int(acTime[0]), int(acTime[1]), int(secsTime[0]))
            self.solarZenith = 90-arcsiUtils.str2Float(headerParams["SunElevationAtCenter"])
            self.solarAzimuth = arcsiUtils.str2Float(headerParams["SunAziumthAtCenter"])
            # Get the geographic lat/long corners of the image.
            self.latTL = arcsiUtils.str2Float(headerParams["ProdULLat"])
            self.lonTL = arcsiUtils.str2Float(headerParams["ProdULLon"])
            self.latTR = arcsiUtils.str2Float(headerParams["ProdURLat"])
            self.lonTR = arcsiUtils.str2Float(headerParams["ProdURLon"])
            self.latBL = arcsiUtils.str2Float(headerParams["ProdLLLat"])
            self.lonBL = arcsiUtils.str2Float(headerParams["ProdLLLon"])
            self.latBR = arcsiUtils.str2Float(headerParams["ProdLRLat"])
            self.lonBR = arcsiUtils.str2Float(headerParams["ProdLRLon"])
            # Get the projected X/Y corners of the image
            self.xTL = arcsiUtils.str2Float(headerParams["ProdULMapX"])
            self.yTL = arcsiUtils.str2Float(headerParams["ProdULMapY"])
            self.xTR = arcsiUtils.str2Float(headerParams["ProdURMapX"])
            self.yTR = arcsiUtils.str2Float(headerParams["ProdURMapY"])
            self.xBL = arcsiUtils.str2Float(headerParams["ProdLLMapX"])
            self.yBL = arcsiUtils.str2Float(headerParams["ProdLLMapY"])
            self.xBR = arcsiUtils.str2Float(headerParams["ProdLRMapX"])
            self.yBR = arcsiUtils.str2Float(headerParams["ProdLRMapY"])
            
            
          
            #Get projection
            inProj = osr.SpatialReference()
            if (headerParams["MapProjection"] == "UTM") and (headerParams["Datum"] == "WGS84") and (headerParams["Ellipsoid"] == "WGS_84"):
                utmZone = int(headerParams["ZoneNo"])
                utmCode = "WGS84UTM" + str(utmZone) + str("N")
                #print("UTM: ", utmCode)
                inProj.ImportFromEPSG(self.epsgCodes[utmCode])
            elif (headerParams["MapProjection"] == "PS") and (headerParams["Datum"] == "WGS84") and (headerParams["Ellipsoid"] == "WGS_84"):
                inProj.ImportFromWkt("PROJCS[\"PS WGS84\", GEOGCS[\"WGS 84\",Datum[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563, AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Polar_Stereographic\"],PARAMETER[\"latitude_of_origin\",-71],PARAMETER[\"central_meridian\",0],PARAMETER[\"scale_factor\",1],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]]]")
            else:
                raise ARCSIException("Expecting LISS III to be projected in UTM or PolarStereographic (PS) with datum=WGS84 and ellipsoid=WGS84.")
            if self.inWKT is "":
                self.inWKT = inProj.ExportToWkt()
                # Check image is square!
            if not ((self.xTL == self.xBL) and (self.yTL == self.yTR) and (self.xTR == self.xBR) and (self.yBL == self.yBR)):
                raise ARCSIException("Image is not square in projected coordinates.")
            self.xCentre = self.xTL + ((self.xTR - self.xTL)/2)
            self.yCentre = self.yBR + ((self.yTL - self.yBR)/2)
            wgs84latlonProj = osr.SpatialReference()
            wgs84latlonProj.ImportFromEPSG(4326)
            wktPt = 'POINT(%s %s)' % (self.xCentre, self.yCentre)
            #print(wktPt)
            point = ogr.CreateGeometryFromWkt(wktPt)
            point.AssignSpatialReference(inProj)
            point.TransformTo(wgs84latlonProj)
            #print(point)
            self.latCentre = point.GetY()
            self.lonCentre = point.GetX()
            #print("Lat: " + str(self.latCentre) + " Long: " + str(self.lonCentre))
            filesDIR = os.path.dirname(inputHeader)
            #getting bandfiles
            self.band3File = os.path.join(filesDIR,"BAND2.tif")
            self.band4File = os.path.join(filesDIR,"BAND3.tif" )
            self.band5File = os.path.join(filesDIR,"BAND4.tif")
            self.band6File = os.path.join(filesDIR,"BAND5.tif")
            
            
           
            self.b3RadMulti = arcsiUtils.str2Float(headerParams["B2SaturationRadiance"])/4095
            self.b4RadMulti = arcsiUtils.str2Float(headerParams["B3SaturationRadiance"])/4095
            self.b5RadMulti = arcsiUtils.str2Float(headerParams["B4SaturationRadiance"])/4095
            self.b6RadMulti = arcsiUtils.str2Float(headerParams["B5SaturationRadiance"])/4095
            
            
            self.b3RadAdd = 0.00
            self.b4RadAdd = 0.00
            self.b5RadAdd = 0.00
            self.b6RadAdd = 0.00
           
            #getting max and min quatizations
           
            self.b3CalMin = 1.00
            self.b3CalMax = 1023.00
            self.b4CalMin = 1.00
            self.b4CalMax = 1023.00
            self.b5CalMin = 1.00
            self.b5CalMax = 1023.00
            self.b6CalMin = 1.00
            self.b6CalMax = 1023.00
            #getting max and min radiance
            self.b3MinRad = arcsiUtils.str2Float(headerParams["B2_Lmin"])*10
            self.b3MaxRad = arcsiUtils.str2Float(headerParams["B2_Lmax"])*10
            self.b4MinRad = arcsiUtils.str2Float(headerParams["B3_Lmin"])*10
            self.b4MaxRad = arcsiUtils.str2Float(headerParams["B3_Lmax"])*10
            self.b5MinRad = arcsiUtils.str2Float(headerParams["B4_Lmin"])*10
            self.b5MaxRad = arcsiUtils.str2Float(headerParams["B4_Lmax"])*10
            self.b6MinRad = arcsiUtils.str2Float(headerParams["B5_Lmin"])*10
            self.b6MaxRad = arcsiUtils.str2Float(headerParams["B5_Lmax"])*10
            #calculating earth sun disatnce
            dist = datetime.datetime.strptime(headerParams["SceneCenterTime"], '%d-%b-%Y %H:%M:%S.%f')
            self.earthSunDistance = coordinates.sun.earth_distance(dist)
            #getting resolution   
            self.gridCellSizeRefl = arcsiUtils.str2Float(headerParams["OutputResolutionAlong"], 24.0)
            #getting file generation date and time
            fileDateStr = headerParams["GenerationDateTime"].strip()
            self.fileDateObj = datetime.datetime.strptime(fileDateStr, "%d-%b-%Y %H:%M:%S")
            self.satAltitude = arcsiUtils.str2Float(headerParams["SatelliteAltitude"])
            

        except Exception as e:
            raise e

    def getSolarIrrStdSolarGeom(self):
        """
        Get Solar Azimuth and Zenith as standard geometry.
        Azimuth: N=0, E=90, S=180, W=270.
        """
        solarAz = rsgislib.imagecalibration.solarangles.getSolarIrrConventionSolarAzimuthFromUSGS(self.solarAzimuth)
        return (solarAz, self.solarZenith)
    

    def getSensorViewGeom(self):
        """
        Get sensor viewing angles
        returns (viewAzimuth, viewZenith)
        """
        return (0.0, 0.0)

    def generateOutputBaseName(self):
        """
        Provides an implementation for the LISS III sensor
        """
        rowpath = "r" + str(self.row) + "p" + str(self.path)
        outname = self.defaultGenBaseOutFileName()
        outname = outname + str("_") + rowpath
        return outname

    def generateMetaDataFile(self, outputPath, outputFileName, productsStr, validMaskImage="", footprintCalc=False, calcdValuesDict=dict(), outFilesDict=dict()):
        """
        Generate file metadata.
        """
        outJSONFilePath = os.path.join(outputPath, outputFileName)
        jsonData = self.getJSONDictDefaultMetaData(productsStr, validMaskImage, footprintCalc, calcdValuesDict, outFilesDict)
        sensorInfo = jsonData['SensorInfo']
        sensorInfo['Row'] = self.row
        sensorInfo['Path'] = self.path
        sensorInfo['SensorID'] = self.sensorID
        sensorInfo['SpacecraftID'] = self.spacecraftID
        acqDict = jsonData['AcquasitionInfo']
        #acqDict['EarthSunDistance'] = self.earthSunDistance
        imgInfo = dict()                     
        imgInfo['CellSizeRefl'] = self.gridCellSizeRefl        
        jsonData['ImageInfo'] = imgInfo
        with open(outJSONFilePath, 'w') as outfile:
            json.dump(jsonData, outfile, sort_keys=True,indent=4, separators=(',', ': '), ensure_ascii=False)

    def expectedImageDataPresent(self):
        imageDataPresent = True       
        if not os.path.exists(self.band3File):
            imageDataPresent = False
        if not os.path.exists(self.band4File):
            imageDataPresent = False
        if not os.path.exists(self.band5File):
            imageDataPresent = False
        if not os.path.exists(self.band6File):
            imageDataPresent = False      
        return imageDataPresent

    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat, outWKTFile):
        raise ARCSIException("LISS III  does not provide any image masks, do not use the MASK option.")

    def mosaicImageTiles(self, outputPath):
        raise ARCSIException("Image data does not need mosaicking")

    def resampleImgRes(self, outputPath, resampleToLowResImg, resampleMethod='cubic', multicore=False):
        raise ARCSIException("Image data does not need resampling")

    def sharpenLowResRadImgBands(self, inputImg, outputImage, outFormat):
        raise ARCSIException("Image sharpening is not available for this sensor.")

    def generateValidImageDataMask(self, outputPath, outputMaskName, viewAngleImg, outFormat):
        print("Create the valid data mask")
        tmpBaseName = os.path.splitext(outputMaskName)[0]
        tmpValidPxlMsk = os.path.join(outputPath, tmpBaseName+'vldpxlmsk.kea')
        outputImage = os.path.join(outputPath, outputMaskName)
        inImages = [self.band3File, self.band4File, self.band5File,self.band6File]
        rsgislib.imageutils.genValidMask(inimages=inImages, outimage=tmpValidPxlMsk, gdalformat='KEA', nodata=0.0)
        rsgislib.rastergis.populateStats(tmpValidPxlMsk, True, False, True)
        # Check there is valid data
        ratDS = gdal.Open(tmpValidPxlMsk, gdal.GA_ReadOnly)
        Histogram = rat.readColumn(ratDS, "Histogram")
        ratDS = None
        if Histogram.shape[0] < 2:
            raise ARCSIException("There is no valid data in this image.")
        if not os.path.exists(viewAngleImg):
            print("Calculate Image Angles")
            rsgislib.rastergis.spatialExtent(clumps=tmpValidPxlMsk, minXX='MinXX', minXY='MinXY', maxXX='MaxXX', maxXY='MaxXY', minYX='MinYX', minYY='MinYY', maxYX='MaxYX', maxYY='MaxYY',ratband = 1)
            rsgislib.imagecalibration.calcNadirImgViewAngle(tmpValidPxlMsk, viewAngleImg, 'KEA', 821365, 'MinXX', 'MinXY', 'MaxXX', 'MaxXY', 'MinYX', 'MinYY', 'MaxYX', 'MaxYY')
        rsgislib.imagecalc.imageMath(viewAngleImg, outputImage, 'b1<26?1:0', outFormat, rsgislib.TYPE_64FLOAT)
        rsgisUtils = rsgislib.RSGISPyUtils()
        rsgisUtils.deleteFileWithBasename(tmpValidPxlMsk)
        return outputImage

    def convertImageToRadiance(self, outputPath, outputReflName, outputThermalName, outFormat):
        print("Converting to Radiance")
        outputImage = os.path.join(outputPath, outputReflName)
        bandDefnSeq = list()
        lsBand = collections.namedtuple('LSBand', ['bandName', 'fileName', 'bandIndex', 'lMin', 'lMax', 'qCalMin', 'qCalMax'])
        bandDefnSeq.append(lsBand(bandName="Green", fileName=self.band3File, bandIndex=1, lMin=self.b3MinRad, lMax=self.b3MaxRad, qCalMin=self.b3CalMin, qCalMax=self.b3CalMax))
        bandDefnSeq.append(lsBand(bandName="Red", fileName=self.band4File, bandIndex=1, lMin=self.b4MinRad, lMax=self.b4MaxRad, qCalMin=self.b4CalMin, qCalMax=self.b4CalMax))
        bandDefnSeq.append(lsBand(bandName="VNIR", fileName=self.band5File, bandIndex=1, lMin=self.b5MinRad, lMax=self.b5MaxRad, qCalMin=self.b5CalMin, qCalMax=self.b5CalMax))
        bandDefnSeq.append(lsBand(bandName="SWIR", fileName=self.band6File, bandIndex=1, lMin=self.b6MinRad, lMax=self.b6MaxRad, qCalMin=self.b6CalMin, qCalMax=self.b6CalMax))
        rsgislib.imagecalibration.landsat2Radiance(outputImage, outFormat, bandDefnSeq)
        return outputImage, None


    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        lsBand = collections.namedtuple('LSBand', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
              
        bandDefnSeq.append(lsBand(bandName="Green", fileName=self.band3File, bandIndex=1, satVal=self.b3CalMax))
        bandDefnSeq.append(lsBand(bandName="Red", fileName=self.band4File, bandIndex=1, satVal=self.b4CalMax))
        bandDefnSeq.append(lsBand(bandName="VNIR", fileName=self.band5File, bandIndex=1, satVal=self.b5CalMax))
        bandDefnSeq.append(lsBand(bandName="SWIR", fileName=self.band6File, bandIndex=1, satVal=self.b6CalMax))
       
        rsgislib.imagecalibration.saturatedPixelsMask(outputImage, outFormat, bandDefnSeq)

        return outputImage


    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        solarIrradianceVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        solarIrradianceVals.append(IrrVal(irradiance=1818.30))
        solarIrradianceVals.append(IrrVal(irradiance=1559.20))
        solarIrradianceVals.append(IrrVal(irradiance=1089.60))
        solarIrradianceVals.append(IrrVal(irradiance=243.80))
               
        rsgislib.imagecalibration.radiance2TOARefl(inputRadImage, outputImage, outFormat, rsgislib.TYPE_64FLOAT, scaleFactor, self.acquisitionTime.year, self.acquisitionTime.month, self.acquisitionTime.day, self.solarZenith, solarIrradianceVals)
        return outputImage
    
    def defineDarkShadowImageBand(self):
        return 3       
    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        sixsCoeffs = numpy.zeros((4, 6), dtype=numpy.float32)
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        #s.ground_reflectance = Py6S.GroundReflectance.HomogeneousHapke(0.101, -0.263, 0.589, 0.046)
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.User()
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute)/60.0
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(surfaceAltitude)
        s.altitudes.set_sensor_satellite_level()
        if useBRDF:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrBRDFFromRadiance(200)
        else:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal
       
        # Band 3
        s.wavelength = Py6S.Wavelength(0.52, 0.59, [ 0.043 ,  0.1695,  0.385 ,  0.566 ,  0.626 ,  0.6225,  0.644 ,
                0.6895,  0.726 ,  0.7162,  0.7005,  0.711 ,  0.732 ,  0.736 ,
                0.749 ,  0.7792,  0.797 ,  0.7792,  0.7885,  0.857 ,  0.941 ,
                0.9327,  0.9135,  0.9377,  0.991 ,  0.9948,  0.9615,  0.8802,
                0.721 ])
        s.run()
        sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 4
        s.wavelength = Py6S.Wavelength(0.62, 0.68, [0.4751,  0.7029,  0.8577,  0.9173,  0.917 ,  0.8958,  0.8786,
                    0.8765,  0.8929,  0.9241,  0.957 ,  0.98  ,  0.9854,  0.9671,
                    0.9421,  0.9251,  0.9229,  0.9375,  0.9623,  0.9857,  0.9981,
                    0.9893,  0.9731,  0.9813,  0.9933])
        s.run()
        sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 5
        s.wavelength = Py6S.Wavelength(0.77, 0.86, [0.595 ,  0.6905,  0.825 ,  0.916 ,  0.974 ,  0.9955,  0.9945,
                    0.9815,  0.967 ,  0.946 ,  0.9185,  0.891 ,  0.871 ,  0.845 ,
                    0.8285,  0.8113,  0.796 ,  0.7805,  0.768 ,  0.7562,  0.747 ,
                    0.7342,  0.7325,  0.727 ,  0.723 ,  0.7162,  0.7095,  0.702 ,
                    0.692 ,  0.677 ,  0.6585,  0.628 ,  0.599 ,  0.5645,  0.5175,
                    0.4698,  0.413])
        s.run()
        sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 6
        s.wavelength = Py6S.Wavelength(1.55, 1.70, [0.404 ,  0.4703,  0.529 ,  0.5982,  0.667 ,  0.6973,  0.751 ,
                    0.798 ,  0.837 ,  0.87  ,  0.898 ,  0.9208,  0.938 ,  0.9515,
                    0.9625,  0.9712,  0.977 ,  0.9828,  0.986 ,  0.9885,  0.99  ,
                    0.9925,  0.9785,  0.954 ,  0.959 ,  0.949 ,  0.9195,  0.9028,
                    0.887 ,  0.8738,  0.8705,  0.8387,  0.816 ,  0.8023,  0.7885,
                    0.777 ,  0.766 ,  0.7565,  0.749 ,  0.7408,  0.735 ,  0.73  ,
                    0.727 ,  0.7248,  0.721 ,  0.7248,  0.7065,  0.692 ,  0.687 ,
                    0.6728,  0.6475,  0.6132,  0.586 ,  0.5558,  0.5235,  0.4882,
                    0.449 ,  0.381 ,  0.363 ,  0.3172,  0.272])
        s.run()
        sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[3,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[3,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[3,5] = float(s.outputs.values['environmental_irradiance'])
      
        return sixsCoeffs

    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF, scaleFactor):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)
        imgBandCoeffs = list()
        sixsCoeffs = self.calc6SCoefficients(aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF)        
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
        imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))

        rsgislib.imagecalibration.apply6SCoeffSingleParam(inputRadImage, outputImage, outFormat, rsgislib.TYPE_64FLOAT, scaleFactor, 0, True, imgBandCoeffs)
        return outputImage

    def convertImageToSurfaceReflDEMElevLUT(self, inputRadImage, inputDEMFile, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, scaleFactor, elevCoeffs=None):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevCoeffs is None:
            print("Build an LUT for elevation values.")
            elev6SCoeffsLUT = self.buildElevation6SCoeffLUT(aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax)
            print("LUT has been built.")
            elevCoeffs = list()
            for elevLUT in elev6SCoeffsLUT:
                imgBandCoeffs = list()
                sixsCoeffs = elevLUT.Coeffs
                elevVal = elevLUT.Elev                
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
                imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))
               
                elevCoeffs.append(rsgislib.imagecalibration.ElevLUTFeat(Elev=float(elevVal), Coeffs=imgBandCoeffs))

        rsgislib.imagecalibration.apply6SCoeffElevLUTParam(inputRadImage, inputDEMFile, outputImage, outFormat, rsgislib.TYPE_64FLOAT, scaleFactor, 0, True, elevCoeffs)
        return outputImage, elevCoeffs


    def convertImageToSurfaceReflAOTDEMElevLUT(self, inputRadImage, inputDEMFile, inputAOTImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax, scaleFactor, elevAOTCoeffs=None):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevAOTCoeffs is None:
            print("Build an LUT for elevation and AOT values.")
            elevAOT6SCoeffsLUT = self.buildElevationAOT6SCoeffLUT(aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax)
            elevAOTCoeffs = list()
            for elevLUT in elevAOT6SCoeffsLUT:
                elevVal = elevLUT.Elev
                aotLUT = elevLUT.Coeffs
                aot6SCoeffsOut = list()
                for aotFeat in aotLUT:
                    sixsCoeffs = aotFeat.Coeffs
                    aotVal = aotFeat.AOT
                    imgBandCoeffs = list()
                   
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
                    imgBandCoeffs.append(rsgislib.imagecalibration.Band6SCoeff(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))               
                    aot6SCoeffsOut.append(rsgislib.imagecalibration.AOTLUTFeat(AOT=float(aotVal), Coeffs=imgBandCoeffs))
                elevAOTCoeffs.append(rsgislib.imagecalibration.ElevLUTFeat(Elev=float(elevVal), Coeffs=aot6SCoeffsOut))
        rsgislib.imagecalibration.apply6SCoeffElevAOTLUTParam(inputRadImage, inputDEMFile, inputAOTImage, outputImage, outFormat, rsgislib.TYPE_64FLOAT, scaleFactor, 0, True, elevAOTCoeffs)

        return outputImage, elevAOTCoeffs

    def run6SToOptimiseAODValue(self, aotVal, radBlueVal, predBlueVal, aeroProfile, atmosProfile, grdRefl, surfaceAltitude):
        """Used as part of the optimastion for identifying values of AOD"""
        print("Testing AOD Val: ", aotVal,)
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.User()
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute)/60.0
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(surfaceAltitude)
        s.altitudes.set_sensor_satellite_level()
        s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal
        # Band 2 (Blue!)
        s.wavelength = Py6S.Wavelength(0.436, 0.5285, [0.000010, 0.000117, 0.000455, 0.001197, 0.006869, 0.027170, 0.271370, 0.723971, 0.903034, 0.909880, 0.889667, 0.877453, 0.879688, 0.891913, 0.848533, 0.828339, 0.868497, 0.912538, 0.931726, 0.954248, 0.956424, 0.978564, 0.989469, 0.968801, 0.988729, 0.967361, 0.966125, 0.981834, 0.963135, 0.996498, 0.844893, 0.190738, 0.005328, 0.001557, 0.000516, 0.000162, 0.000023, -0.000016])
        s.run()
        aX = float(s.outputs.values['coef_xa'])
        bX = float(s.outputs.values['coef_xb'])
        cX = float(s.outputs.values['coef_xc'])
        tmpVal = (aX*radBlueVal)-bX;
        reflBlueVal = tmpVal/(1.0+cX*tmpVal)
        outDist = math.sqrt(math.pow((reflBlueVal - predBlueVal),1))
        print("\taX: ", aX, " bX: ", bX, " cX: ", cX, "     Dist = ", outDist)
        return outDist
   
    def estimateImageToAODUsingDOS(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax, globalDOS, simpleDOS, dosOutRefl):
        try:
            print("Estimating AOD Using DOS")
            arcsiUtils = ARCSIUtils()
            outputAOTImage = os.path.join(outputPath, outputName)
            tmpBaseName = os.path.splitext(outputName)[0]
            imgExtension = arcsiUtils.getFileExtension(outFormat)
            dosBlueImage = ""
            minObjSize = 3
            darkPxlPercentile = 0.01
            blockSize = 1000
            if simpleDOS:
                outputDOSBlueName = tmpBaseName + "DOSBlue" + imgExtension
                dosBlueImage, bandOff = self.convertImageBandToReflectanceSimpleDarkSubtract(inputTOAImage, outputPath, outputDOSBlueName, outFormat, dosOutRefl, 1)
            elif globalDOS:
                dosBlueImage = self.performDOSOnSingleBand(inputTOAImage, 1, outputPath, tmpBaseName, "Blue", "KEA", tmpPath, minObjSize, darkPxlPercentile, dosOutRefl)
            else:
                dosBlueImage = self.performLocalDOSOnSingleBand(inputTOAImage, 1, outputPath, tmpBaseName, "Blue", "KEA", tmpPath, minObjSize, darkPxlPercentile, blockSize, dosOutRefl)

            thresImageClumpsFinal = os.path.join(tmpPath, tmpBaseName + "_clumps" + imgExtension)
            rsgislib.segmentation.segutils.runShepherdSegmentation(inputTOAImage, thresImageClumpsFinal, tmpath=tmpPath, gdalformat="KEA", numClusters=40, minPxls=10, bands=[1,2,3,4], processInMem=True)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanElev"))
            rsgislib.rastergis.populateRATWithStats(inputDEMFile, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcTOA = list()
            stats2CalcTOA.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanB2DOS"))
            rsgislib.rastergis.populateRATWithStats(dosBlueImage, thresImageClumpsFinal, stats2CalcTOA)

            stats2CalcRad = list()
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=1, meanField="MeanB2RAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=3, meanField="MeanB4RAD"))
            stats2CalcRad.append(rsgislib.rastergis.BandAttStats(band=2, meanField="MeanB3RAD"))
            
            rsgislib.rastergis.populateRATWithStats(inputRADImage, thresImageClumpsFinal, stats2CalcRad)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            Histogram = rat.readColumn(ratDS, "Histogram")
            MeanElev = rat.readColumn(ratDS, "MeanElev")
            
#           
            
            
            MeanB3RAD = rat.readColumn(ratDS, "MeanB3RAD")
            MeanB4RAD = rat.readColumn(ratDS, "MeanB4RAD")            
            radNDVI = (MeanB4RAD - MeanB3RAD)/(MeanB4RAD + MeanB3RAD)
            
            selected = Histogram * 2
            selected[...] = 0
            selected[radNDVI>0.2] = 1
            rat.writeColumn(ratDS, "Selected", selected)
            ratDS = None

            rsgislib.rastergis.spatialLocation(thresImageClumpsFinal, "Eastings", "Northings")
            rsgislib.rastergis.selectClumpsOnGrid(thresImageClumpsFinal, "Selected", "PredictAOTFor", "Eastings", "Northings", "MeanB2DOS", "min", 20, 20)

            ratDS = gdal.Open(thresImageClumpsFinal, gdal.GA_Update)
            MeanB1DOS = rat.readColumn(ratDS, "MeanB2DOS")
            MeanB1DOS = MeanB1DOS / 1000
            MeanB1RAD = rat.readColumn(ratDS, "MeanB2RAD")
            PredictAOTFor = rat.readColumn(ratDS, "PredictAOTFor")
            numAOTValTests = int(math.ceil((aotValMax - aotValMin)/0.05))+1
            if not numAOTValTests >= 1:
                raise ARCSIException("min and max AOT range are too close together, they need to be at least 0.05 apart.")
            cAOT = aotValMin
            cDist = 0.0
            minAOT = 0.0
            minDist = 0.0
            aotVals = numpy.zeros_like(MeanB1RAD, dtype=numpy.float)
            for i in range(len(MeanB1RAD)):
                if PredictAOTFor[i] == 1:
                    print("Predicting AOD for Segment ", i)
                    for j in range(numAOTValTests):
                        cAOT = aotValMin + (0.05 * j)
                        cDist = self.run6SToOptimiseAODValue(cAOT, MeanB1RAD[i], MeanB1DOS[i], aeroProfile, atmosProfile, grdRefl, MeanElev[i]/1000)
                        if j == 0:
                            minAOT = cAOT
                            minDist = cDist
                        elif cDist < minDist:
                            minAOT = cAOT
                            minDist = cDist
                    #predAOTArgs = (MinB2RAD[i], MeanB2DOS[i], aeroProfile, atmosProfile, grdRefl, MeanElev[i]/1000)
                    #res = minimize(self.run6SToOptimiseAODValue, minAOT, method='nelder-mead', options={'maxiter': 20, 'xtol': 0.001, 'disp': True}, args=predAOTArgs)
                    #aotVals[i] = res.x[0]
                    aotVals[i] = minAOT
                    print("IDENTIFIED AOT: ", aotVals[i])
                else:
                    aotVals[i] = 0
            rat.writeColumn(ratDS, "AOT", aotVals)

            Eastings = rat.readColumn(ratDS, "Eastings")
            Northings = rat.readColumn(ratDS, "Northings")
            ratDS = None

            Eastings = Eastings[PredictAOTFor!=0]
            Northings = Northings[PredictAOTFor!=0]
            aotVals = aotVals[PredictAOTFor!=0]

            interpSmoothing = 10.0
            self.interpolateImageFromPointData(inputTOAImage, Eastings, Northings, aotVals, outputAOTImage, outFormat, interpSmoothing, True, 0.05)
            if not self.debugMode:
                gdalDriver = gdal.GetDriverByName(outFormat)
                gdalDriver.Delete(thresImageClumpsFinal)
                gdalDriver.Delete(dosBlueImage)
            return outputAOTImage
        except Exception as e:
            raise e

    def estimateSingleAOTFromDOS(self, radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl):
        try:
            return self.estimateSingleAOTFromDOSBandImpl(radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl, 1)
        except Exception as e:
            raise

    def setBandNames(self, imageFile):
        dataset = gdal.Open(imageFile, gdal.GA_Update)       
        dataset.GetRasterBand(1).SetDescription("Green")
        dataset.GetRasterBand(2).SetDescription("Red")
        dataset.GetRasterBand(3).SetDescription("VNIR")
        dataset.GetRasterBand(4).SetDescription("SWIR")        
        dataset = None

    def cleanLocalFollowProcessing(self):
        print("")



