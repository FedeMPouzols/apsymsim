# THE SYMBOL "#" IS USED TO INSERT COMMENTS, 
# WHICH WILL NOT BE READ BY THE PROGRAM


# THIS IS THE IMAGE SIZE (IN ARCSEC.):

IMSIZE 3.0


# THIS IS THE OBSERVING WAVELENGTH (IN METERS):

WAVELENGTH  1.e-03


# HERE COME THE MODEL COMPONENTS. THERE ARE THREE
# RECOGNIZED KINDS OF COMPONENTS:
# POINT SOURCE (P)
# GAUSSIAN     (G)
# UNIFORM DISC (D)
#
# EACH COMPONENT IS SET IN ONE LINE WITH 5 ROWS:
#
# KIND, RA OFFSET, DEC OFFSET, FLUX DENSITY, SIZE
#
# WHERE "KIND" IS EITHER P, G, OR D;
# RA AND DEC OFFSETS ARE IN ARCSEC.
# FLUX DENSITY IS IN JY.
# SIZE IS IN ARCSEC (ONLY GIVEN FOR DISCS AND GAUSSIANS).

P  -0.3  0.0  1.0
P   0.3  0.0  1.0

