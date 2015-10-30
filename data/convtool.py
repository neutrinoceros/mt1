#-*-coding:utf-8-*-
GSI =  6.67408 * 1e-11# m^3 kg^-1 s^-2

m2au    = 1.49597870700e11 #from the same article we get our data from
kg2Msol = 132712440041.939400 * 1e9 /GSI #idem
s2day   = 3600 * 24
si2sad  = m2au**3 * kg2Msol**-1 * s2day**-2

GSAD =  GSI * si2sad

#print si2sad, GSAD

# esi = GSI*1.*1./1. #energie pot en si de deux objets de 1kg séparés par 1m 
# esad = GSAD*(kg2Msol**-1 * 1.)**2/(m2au**-1) #la même en sad
# print esi, esad
#SAD stands for "Solar mass, Astronomical unit, Day"
