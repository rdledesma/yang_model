# -*- coding: utf-8 -*-
"""
Created on Mon May  31 17:08:29 2022

@author: Dario Ledesma
"""
import math
import pandas as pd
import numpy as np


class Geo:
    
    def __init__(self, range_dates,  lat, long, gmt, alt, beta ):
        self.gmt = gmt
        self.lat = lat
        self.long = long
        self.altura = alt
        self.df = pd.DataFrame()
        
        """
        genera las filas para todo el df
        dado fechade de inicio y fin espaciado por freq cant de tiempo
        """
        #self.df['Fecha'] =  pd.date_range(start=desde+' 00:00:00', end=hasta+' 23:59:00', freq= freq+' min')
        self.df['Fecha'] = range_dates
        """
        #el desplamiento para realizar 
        #el cálculo en el centro del invervalo
        """
        #self.df['Fecha'] =  self.df['Fecha']
        self.df['N'] = self.df['Fecha'].dt.day_of_year
        self.df['M'] = self.df['Fecha'].dt.month
        self.df['Y'] = self.df['Fecha'].dt.year
        self.df['E'] = self.df.apply(lambda r: self.getE(r['N'], r['Y']), axis=1)
        self.df['HR'] = self.df['Fecha'].dt.hour + (self.df['Fecha'].dt.minute) /60 + (self.df['Fecha'].dt.second)/3600
        self.df['HS'] = self.getHS()

        self.df['delta rad'] = self.df.apply(lambda r: self.delta(r['N'], r['Y']), axis=1)
        self.df['delta'] = self.df['delta rad'].apply(math.degrees)
        self.df['ws'] = self.df['delta rad'].apply(self.getWs)
        
        self.df['Fn'] = self.df.apply(lambda r: self.Fn(r['N'], r['Y']), axis=1)
        
          
        self.df['w'] = 15 * (12 - self.df['HS'])    
        self.df['w rad'] = self.df['w'].apply(math.radians)

        self.df['CTZ'] =  self.df.apply(lambda r: self.getCTZ(r['delta rad'], self.lat, r['w rad']), axis=1)

        self.df['TZ'] = self.df['CTZ'].apply(math.acos)
 
        
        self.df['Ys'] = self.df.apply(lambda r: self.Ys(r['w'], r['CTZ'], r['TZ'], r['delta rad']), axis=1)
        
        self.df['CT'] = np.where(self.df.CTZ>0, self.df.CTZ * math.cos(math.radians(beta)) + self.df.TZ.apply(math.sin) * math.sin(math.radians(beta)) * self.df.Ys.apply(math.cos), 0)
        
        self.df['Ys'] = self.df.apply(lambda r: self.Ys(r['w'], r['CTZ'], r['TZ'], r['delta rad']), axis=1)
        self.df['alphaS'] = self.df['CTZ'].apply(math.asin).apply(math.degrees)
        self.df['E0'] = self.df.apply(lambda r: self.getE0(r['N'], r['Y']), axis=1)
        self.df['TOA'] = self.df.apply(lambda r: self.TOA(r['E0'], r['CTZ']), axis=1 )
        self.df['I0'] = self.df.apply(lambda r: 
                                      self.TOADiaria(
                                          r['E0'],
                                          r['delta'],
                                          r['CTZ'],
                                          r['ws']
                                                      ), axis=1)
        #self.df['Ma'] = self.generateMa()
        #self.df['Ma2'] = self.getMA(self.df['CTZ'])
        
        ##ARGPV1
        ##Masa de aire de casten, es la que se utiliza
        self.df['Mak'] = self.MaK()
        self.ktrp = 0.7 + 1.6391 * 10**-3 * self.altura ** 0.5500 
        
        #self.ktrp_2 = 0.649 + 0.02 * self.altura ** 0.28
        
        self.df['GHIargp'] = self.generateGHIargp(self.df)
        # self.df['GHIargp_2'] = self.generateGHIargp_2(self.df)
        
        
        
    #Ecuación de tiempo
    #dado dia ordinal 
    #devuelve valor ecuación del tiempo en minutos
    def getE(self,n,y):
        
        if(y%4) == 0:
            gamma = 2 * math.pi * (n-1)/366
        else:
            gamma = 2 * math.pi * (n-1)/365
            
        E = 229.18 * (0.000075+ 0.001868 * math.cos(gamma) - 0.032077*math.sin(gamma) - 0.014615*math.cos(2*gamma) - 0.04089*math.sin(2*gamma))
        return E

    def Fn(self, n, y):
            
        if(y%4) == 0:
            gamma = 2 * math.pi * (n-1)/366
        else:
            gamma = 2 * math.pi * (n-1)/365
            
        fn = 1.000110 + 0.034221*math.cos(gamma) + 0.001280*math.sin(gamma) + 0.000719*math.cos(2*gamma) + 0.000077*math.sin(2*gamma)
        return fn
    
    #Hora solar 
    #devuelve valor hora solar
    #recorre cada fila del df
    def getHS(self):
        A = 1
        if self.gmt<0:
            A = -1
        return self.df['HR'] + (4 * ((A * 15 * self.gmt)- (A*self.long))+ self.df['E'])/60
    
  
    
    #Declinación solar
    #dado dia ordinal
    #devuelve declinacion en radianes
    def delta(self, n, y):
        
        if(y%4) == 0:
            gamma = 2 * math.pi * (n-1)/366
        else:
            gamma = 2 * math.pi * (n-1)/365
            
        delta = 0.006918 - 0.399912 * math.cos(gamma) + 0.070257 * math.sin(gamma) - 0.006758 * math.cos(2*gamma) + 0.000907*math.sin(2*gamma) - 0.002697*math.cos(3*gamma)+ 0.00148*math.sin(3*gamma)
        return delta
    
    #CTZ
    #dado decliancio, lat y angulo horario
    #devuelve cos tita z en radianes
    def getCTZ(self, delta, lat, omega):
        latR = math.radians(self.lat)    
        return (math.cos(latR) * math.cos(delta)* math.cos(omega)) + (math.sin(latR)*math.sin(delta))
        #return math.sin(delta) * math.sin(math.radians(lat)) + math.cos(delta) * math.cos(math.radians(lat))* math.cos(omega)


    def getWs(self, delta):
        return math.acos(-math.tan(delta) * (math.tan(math.radians(self.lat))))


    #E0
    #dado dia ordinal
    #devuelve factor de correción orbital
    def getE0(self, N, y):
        if(y%4) == 0:
            return 1+0.033* math.cos(2* math.pi * N /366)
        else:
            return 1+0.033* math.cos(2* math.pi * N /365)
            
        
        
    
    #T0A
    #dado fac. corr. orb. y ctz
    #devuelve irrad. solar extr. expr en whm2
    def TOA(self, E0, CTZ):
        if CTZ<0:
            return 0
        else:
            return 1361 * E0 * CTZ
        
        
    ##Ys
    def Ys(self, omega, CTZ, TZ, d):
        
        
        sinTitaZ = math.sin(TZ)
        sinLat = math.sin(math.radians(self.lat))
        cosLat = math.cos(math.radians(self.lat))
        signow = 1
        
        if(omega<0):
            signow = -1
        
        
        #ys = math.acos((math.sin(d) - CTZ * math.sin(math.radians(self.lat))) / (sinTitaZ * math.cos(math.radians(self.lat))))
        
        ys = math.acos((CTZ * sinLat - math.sin(d)) / (sinTitaZ * cosLat ) )
        
        #ACOS((O2*SENO(C2)-SENO(L2))/(SENO(P2)*COS(C2)))
        
        return signow * abs(ys)
    
    #T0A
    #dado fac. corr. orb. y ctz
    #devuelve irrad. solar extr. expr en whm2
    def TOADiaria(self, E0, dlt, CTZ, ws ):
    
        # E0 = data['E0']
        cosDelta = math.cos(dlt)
        sinDelta = math.sin(dlt)
        cosLatR = math.cos(math.radians(self.lat))    
        sinLatR = math.sin(math.radians(self.lat))   
        sinWs = math.sin(ws)
        
        if CTZ>0:
            return abs(24/math.pi *1361*E0 * (cosDelta * cosLatR * sinWs + ws * sinDelta * sinLatR))
        else:
            return 0
        
        
    def getMA(self,CTZs):
        result = []
        for ctz in CTZs:
            try:
                valor1 = 1.002432 * ctz**2 + 0.148386*ctz + 0.0096
                valor2 = ctz**3 + 0.149864*ctz**2 + 0.0102963*ctz + 0.000303978
                result.append(valor1/valor2)
            except Exception:
                result.append(0)
        return result
    
    
    
        
    
    def MaK(self):
        presion = 101355* (288.15/(288.15 - 0.0065 * self.altura)) ** -5.255877
        
        CTZ = self.df['CTZ']
        TZ = self.df['TZ']
        
        Amc = []
        for i, val in enumerate(CTZ):
            Amk = 1/ (val + 0.15*(93.885 - TZ[i])**-1.253)
            Amc.append( Amk * (presion/101355))
            
        return Amc
    
    
    def generateMa(self):
        
        cosTZ = self.df['CTZ'].tolist()
        tz = self.df['TZ'].tolist()
           
        
        
        
        presion =  math.pow(288.15/(288.15 - 0.0065 * 1150), -5.255877);
        results = []
        
        for i, val in enumerate(cosTZ):
            
            try:
                calc = 1/(val + 0.15 * math.pow((93.885-  tz[i]), -1.253));
                results.append(calc * presion)
            except Exception:
                results.append(0)        
        return results
    
   
    
    def generateGHIargp(self, data):
        GHI = data['TOA'].tolist()
        AM = data['Mak'].tolist()
    
    
    
        result = []
        for i, val in enumerate(GHI):
            try:
                result.append(GHI[i] * math.pow( self.ktrp ,math.pow(AM[i], 0.678)))
            except Exception:
                result.append(0)
        return result
    
    def generateGHIargp_2(self, data):
        GHI = data['TOA'].tolist()
        AM = data['Mak'].tolist()
    
    
    
        result = []
        for i, val in enumerate(GHI):
            try:
                result.append(GHI[i] * math.pow( self.ktrp_2 ,math.pow(AM[i], 0.678)))
            except Exception:
                result.append(0)
        return result
    
    
    def to_csv(self, name):
        self.df.to_csv(name+".csv")

 
        
#geo  = Geo(freq='1', lat=-31.28, long=-57.92, gmt=-3, alt=1000, desde='01/01/2016', hasta='01/01/2016').df

