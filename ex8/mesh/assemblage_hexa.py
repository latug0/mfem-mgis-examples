import gmsh
import math
import time
import argparse

def plaque(offset, FuelThickness, CladThickness, RayonCintrage, angle):
    epaissGaineHaut = (PlateHeight - FuelHeight) / 2.
    hautCombB = hautGaineB + epaissGaineHaut     
    hautCombH = hautGaineH - epaissGaineHaut

    rintComb = RayonCintrage - FuelThickness / 2.
    rextComb = RayonCintrage + FuelThickness / 2.
    rintGaine = rintComb - CladThickness
    rextGaine = rextComb + CladThickness
    epBord = StiffenerHalfThickness + CladThickness

    angStifD = angle - AnglePlate / 2.
    angStifG = angle + AnglePlate / 2.

    cosAngStifD = math.cos(2*math.pi*angStifD/360)
    sinAngStifD = math.sin(2*math.pi*angStifD/360)
    tanAngStifD = sinAngStifD/cosAngStifD

    cosAngStifG = math.cos(2*math.pi* angStifG/360)
    sinAngStifG = math.sin(2*math.pi* angStifG/360)
    tanAngStifG = sinAngStifG/cosAngStifG

    # combustible

    racDeltaIntCD  = math.sqrt(cosAngStifD*cosAngStifD*(rintComb*rintComb-epBord*epBord))
    racDeltaExtCD  = math.sqrt(cosAngStifD*cosAngStifD*(rextComb*rextComb-epBord*epBord))
    if cosAngStifD < 0. :
        racDeltaIntCD = - racDeltaIntCD
        racDeltaExtCD = - racDeltaExtCD
    valCorX = epBord*sinAngStifD
    valCorY = epBord*cosAngStifD
    gmsh.model.geo.addPoint(racDeltaIntCD-valCorX, tanAngStifD*racDeltaIntCD+valCorY, hautCombB, meshSizeC, 1+offset)
    gmsh.model.geo.addPoint(racDeltaExtCD-valCorX, tanAngStifD*racDeltaExtCD+valCorY, hautCombB, meshSizeC, 2+offset)
    gmsh.model.geo.addPoint(racDeltaExtCD-valCorX, tanAngStifD*racDeltaExtCD+valCorY, hautCombH, meshSizeC, 3+offset)
    gmsh.model.geo.addPoint(racDeltaIntCD-valCorX, tanAngStifD*racDeltaIntCD+valCorY, hautCombH, meshSizeC, 4+offset)

    racDeltaIntCG = math.sqrt(cosAngStifG*cosAngStifG*(rintComb*rintComb-epBord*epBord))
    racDeltaExtCG = math.sqrt(cosAngStifG*cosAngStifG*(rextComb*rextComb-epBord*epBord))
    if cosAngStifG > 0. :
        racDeltaIntCG = - racDeltaIntCG
        racDeltaExtCG = - racDeltaExtCG
    valCorX = epBord*sinAngStifG
    valCorY = epBord*cosAngStifG
    gmsh.model.geo.addPoint(-racDeltaIntCG+valCorX, -tanAngStifG*racDeltaIntCG-valCorY, hautCombB, meshSizeC, 5+offset)
    gmsh.model.geo.addPoint(-racDeltaExtCG+valCorX, -tanAngStifG*racDeltaExtCG-valCorY, hautCombB, meshSizeC, 6+offset)
    gmsh.model.geo.addPoint(-racDeltaExtCG+valCorX, -tanAngStifG*racDeltaExtCG-valCorY, hautCombH, meshSizeC, 7+offset)
    gmsh.model.geo.addPoint(-racDeltaIntCG+valCorX, -tanAngStifG*racDeltaIntCG-valCorY, hautCombH, meshSizeC, 8+offset)

    gmsh.model.geo.addLine(1+offset, 2+offset, 1+offset)
    gmsh.model.geo.addLine(2+offset, 3+offset, 2+offset)
    gmsh.model.geo.addLine(3+offset, 4+offset, 3+offset)
    gmsh.model.geo.addLine(4+offset, 1+offset, 4+offset)
    gmsh.model.geo.addLine(5+offset, 6+offset, 5+offset)
    gmsh.model.geo.addLine(6+offset, 7+offset, 6+offset)
    gmsh.model.geo.addLine(7+offset, 8+offset, 7+offset)
    gmsh.model.geo.addLine(8+offset, 5+offset, 8+offset)

    gmsh.model.geo.addPoint(0, 0, hautCombB, tag=9+offset) 
    gmsh.model.geo.addPoint(0, 0, hautCombH, tag=10+offset)

    gmsh.model.geo.addCircleArc(1+offset, 9+offset,  5+offset,  9+offset)
    gmsh.model.geo.addCircleArc(2+offset, 9+offset,  6+offset, 10+offset)
    gmsh.model.geo.addCircleArc(3+offset, 10+offset, 7+offset, 11+offset)
    gmsh.model.geo.addCircleArc(4+offset, 10+offset, 8+offset, 12+offset)

    gmsh.model.geo.addCurveLoop([1+offset,  2+offset,  3+offset,   4+offset], 1+offset)
    gmsh.model.geo.addCurveLoop([5+offset,  6+offset,  7+offset,   8+offset], 2+offset)
    gmsh.model.geo.addCurveLoop([1+offset, 10+offset, -5-offset,  -9-offset], 3+offset)
    gmsh.model.geo.addCurveLoop([2+offset, 11+offset, -6-offset, -10-offset], 4+offset)
    gmsh.model.geo.addCurveLoop([3+offset, 12+offset, -7-offset, -11-offset], 5+offset)
    gmsh.model.geo.addCurveLoop([4+offset,  9+offset, -8-offset, -12-offset], 6+offset)

    # gaine

    racDeltaIntD  = math.sqrt(cosAngStifD*cosAngStifD*(rintGaine*rintGaine-StiffenerHalfThickness*StiffenerHalfThickness))
    racDeltaExtD  = math.sqrt(cosAngStifD*cosAngStifD*(rextGaine*rextGaine-StiffenerHalfThickness*StiffenerHalfThickness))
    if cosAngStifD < 0. :
        racDeltaIntD = - racDeltaIntD
        racDeltaExtD = - racDeltaExtD
    valCorX = StiffenerHalfThickness*sinAngStifD
    valCorY = StiffenerHalfThickness*cosAngStifD
    gmsh.model.geo.addPoint(racDeltaIntD-valCorX, tanAngStifD*racDeltaIntD+valCorY, hautGaineB, meshSizeG ,21+offset)
    gmsh.model.geo.addPoint(racDeltaExtD-valCorX, tanAngStifD*racDeltaExtD+valCorY, hautGaineB, meshSizeG ,22+offset)
    gmsh.model.geo.addPoint(racDeltaExtD-valCorX, tanAngStifD*racDeltaExtD+valCorY, hautGaineH, meshSizeG ,23+offset)
    gmsh.model.geo.addPoint(racDeltaIntD-valCorX, tanAngStifD*racDeltaIntD+valCorY, hautGaineH, meshSizeG ,24+offset)

    racDeltaIntG = math.sqrt(cosAngStifG*cosAngStifG*(rintGaine*rintGaine-StiffenerHalfThickness*StiffenerHalfThickness))
    racDeltaExtG = math.sqrt(cosAngStifG*cosAngStifG*(rextGaine*rextGaine-StiffenerHalfThickness*StiffenerHalfThickness))
    if cosAngStifG > 0. :
        racDeltaIntG = - racDeltaIntG
        racDeltaExtG = - racDeltaExtG
    valCorX = StiffenerHalfThickness*sinAngStifG
    valCorY = StiffenerHalfThickness*cosAngStifG
    gmsh.model.geo.addPoint(-racDeltaIntG+valCorX, -tanAngStifG*racDeltaIntG-valCorY, hautGaineB, meshSizeG ,25+offset)
    gmsh.model.geo.addPoint(-racDeltaExtG+valCorX, -tanAngStifG*racDeltaExtG-valCorY, hautGaineB, meshSizeG ,26+offset)
    gmsh.model.geo.addPoint(-racDeltaExtG+valCorX, -tanAngStifG*racDeltaExtG-valCorY, hautGaineH, meshSizeG ,27+offset)
    gmsh.model.geo.addPoint(-racDeltaIntG+valCorX, -tanAngStifG*racDeltaIntG-valCorY, hautGaineH, meshSizeG ,28+offset)

    gmsh.model.geo.addLine(21+offset, 22+offset, 21+offset)
    gmsh.model.geo.addLine(22+offset, 23+offset, 22+offset)
    gmsh.model.geo.addLine(23+offset, 24+offset, 23+offset)
    gmsh.model.geo.addLine(24+offset, 21+offset, 24+offset)
    gmsh.model.geo.addLine(25+offset, 26+offset, 25+offset)
    gmsh.model.geo.addLine(26+offset, 27+offset, 26+offset)
    gmsh.model.geo.addLine(27+offset, 28+offset, 27+offset)
    gmsh.model.geo.addLine(28+offset, 25+offset, 28+offset)

    gmsh.model.geo.addPoint(0, 0, hautGaineB, tag=29+offset)
    gmsh.model.geo.addPoint(0, 0, hautGaineH, tag=30+offset) 

    gmsh.model.geo.addCircleArc(21+offset, 29+offset, 25+offset,29+offset) 
    gmsh.model.geo.addCircleArc(22+offset, 29+offset, 26+offset,30+offset)
    gmsh.model.geo.addCircleArc(23+offset, 30+offset, 27+offset,31+offset)
    gmsh.model.geo.addCircleArc(24+offset, 30+offset, 28+offset,32+offset)

    gmsh.model.geo.addCurveLoop([21+offset, 22+offset,  23+offset,  24+offset],20+offset)
    gmsh.model.geo.addCurveLoop([21+offset, 30+offset, -25-offset, -29-offset],21+offset)
    gmsh.model.geo.addCurveLoop([22+offset, 31+offset, -26-offset, -30-offset],22+offset)
    gmsh.model.geo.addCurveLoop([23+offset, 32+offset, -27-offset, -31-offset],23+offset)
    gmsh.model.geo.addCurveLoop([24+offset, 29+offset, -28-offset, -32-offset],24+offset)

    # pour maillage

    gmsh.model.geo.addLine(1+offset, 21+offset, 33+offset)
    gmsh.model.geo.addLine(2+offset, 22+offset, 34+offset)
    gmsh.model.geo.addLine(3+offset, 23+offset, 35+offset)
    gmsh.model.geo.addLine(4+offset, 24+offset, 36+offset)
    gmsh.model.geo.addLine(5+offset, 25+offset, 37+offset)
    gmsh.model.geo.addLine(6+offset, 26+offset, 38+offset)
    gmsh.model.geo.addLine(7+offset, 27+offset, 39+offset)
    gmsh.model.geo.addLine(8+offset, 28+offset, 40+offset)

    gmsh.model.geo.addCurveLoop([ 21+offset, -34-offset,  -1-offset,  33+offset],25+offset)
    gmsh.model.geo.addCurveLoop([ 34+offset,  22+offset, -35-offset,  -2-offset],26+offset)
    gmsh.model.geo.addCurveLoop([ -3-offset,  35+offset,  23+offset, -36-offset],27+offset)
    gmsh.model.geo.addCurveLoop([-33-offset,  -4-offset,  36+offset,  24+offset],28+offset)
    gmsh.model.geo.addCurveLoop([ 25+offset, -38-offset,  -5-offset,  37+offset],29+offset)
    gmsh.model.geo.addCurveLoop([ 38+offset,  26+offset, -39-offset,  -6-offset],30+offset)
    gmsh.model.geo.addCurveLoop([ -7-offset,  39+offset,  27+offset, -40-offset],31+offset)
    gmsh.model.geo.addCurveLoop([-37-offset,  -8-offset,  40+offset,  28+offset],32+offset)
    gmsh.model.geo.addCurveLoop([ 29+offset,  33+offset,  -9-offset, -37-offset],33+offset)
    gmsh.model.geo.addCurveLoop([ 30+offset,  34+offset, -10-offset, -38-offset],34+offset)
    gmsh.model.geo.addCurveLoop([ 31+offset,  35+offset, -11-offset, -39-offset],35+offset)
    gmsh.model.geo.addCurveLoop([ 32+offset,  36+offset, -12-offset, -40-offset],36+offset)
    gmsh.model.geo.addCurveLoop([ 25+offset,  26+offset,   27+offset,  28+offset],37+offset)

    # volumes

    for su in range(1+offset, 7+offset):
        gmsh.model.geo.addSurfaceFilling([su],su)

    for su in range(20+offset, 38+offset):
        gmsh.model.geo.addSurfaceFilling([su],su)

    gmsh.model.geo.addSurfaceLoop([1+offset,  2+offset,  3+offset,  4+offset,  5+offset,  6+offset] , 1+offset)  
    gmsh.model.geo.addSurfaceLoop([3+offset, 21+offset, 25+offset, 29+offset, 33+offset, 34+offset] ,21+offset)
    gmsh.model.geo.addSurfaceLoop([4+offset, 22+offset, 26+offset, 30+offset, 34+offset, 35+offset] ,22+offset)
    gmsh.model.geo.addSurfaceLoop([5+offset, 23+offset, 27+offset, 31+offset, 35+offset, 36+offset] ,23+offset)
    gmsh.model.geo.addSurfaceLoop([6+offset, 24+offset, 28+offset, 32+offset, 36+offset, 33+offset] ,24+offset)
    gmsh.model.geo.addSurfaceLoop([1+offset, 20+offset, 25+offset, 26+offset, 27+offset, 28+offset] ,25+offset)
    gmsh.model.geo.addSurfaceLoop([2+offset, 37+offset, 29+offset, 32+offset, 31+offset, 30+offset] ,26+offset)

    gmsh.model.geo.addVolume([1+offset],1+offset) #comb

    for vo in range(21+offset, 27+offset): #gaine    
        gmsh.model.geo.addVolume([vo],vo)

    gmsh.model.geo.synchronize()


def stiffener(offset, rintStif, rextStif, hautStifB, hautStifH, angle, i, j):
    first = False
    if (i == 0):
        first = True
    
    cosAngStif = math.cos(2*math.pi*angle/360)
    sinAngStif = math.sin(2*math.pi*angle/360)

    valCorX = StiffenerHalfThickness*sinAngStif
    valCorY = StiffenerHalfThickness*cosAngStif
    if first:
        gmsh.model.geo.addPoint(rintStif*cosAngStif+valCorX, rintStif*sinAngStif-valCorY, hautStifB, meshSizeS, 41+offset)
    gmsh.model.geo.addPoint(rextStif*cosAngStif+valCorX, rextStif*sinAngStif-valCorY, hautStifB, meshSizeS, 42+offset)
    gmsh.model.geo.addPoint(rextStif*cosAngStif+valCorX, rextStif*sinAngStif-valCorY, hautStifH, meshSizeS, 43+offset)
    if first:
        gmsh.model.geo.addPoint(rintStif*cosAngStif+valCorX, rintStif*sinAngStif-valCorY, hautStifH, meshSizeS, 44+offset)
        gmsh.model.geo.addPoint(rintStif*cosAngStif-valCorX, rintStif*sinAngStif+valCorY, hautStifB, meshSizeS, 45+offset)
    gmsh.model.geo.addPoint(rextStif*cosAngStif-valCorX, rextStif*sinAngStif+valCorY, hautStifB, meshSizeS, 46+offset)
    gmsh.model.geo.addPoint(rextStif*cosAngStif-valCorX, rextStif*sinAngStif+valCorY, hautStifH, meshSizeS, 47+offset)
    if first:
        gmsh.model.geo.addPoint(rintStif*cosAngStif-valCorX, rintStif*sinAngStif+valCorY, hautStifH, meshSizeS, 48+offset)

    # cas ou ce n'est pas la 1ere couronne, on re-utilise la surface du stiffener
    if first:
        offset2 = offset
    else :
        offset2 = offset - PlatesNbr * offsetBase

    # pour la surface qui fait la jonction stiffener-gaine 
    offsetl = (i * PlatesNbr + (j + 1)%PlatesNbr) * offsetBase

    if first:
        gmsh.model.geo.addLine(41+offset2, 42+offset, 41+offset)
        gmsh.model.geo.addLine(43+offset, 44+offset2, 43+offset)
        gmsh.model.geo.addLine(45+offset2, 46+offset, 45+offset)
        gmsh.model.geo.addLine(47+offset, 48+offset2, 47+offset)
        gmsh.model.geo.addLine(21+offsetl, 45+offset2, 53+offset)
        gmsh.model.geo.addLine(24+offsetl, 48+offset2, 56+offset)
        gmsh.model.geo.addLine(25+offset, 41+offset2, 61+offset)
        gmsh.model.geo.addLine(28+offset, 44+offset2, 64+offset)

        gmsh.model.geo.addLine(44+offset, 41+offset, 44+offset)
        gmsh.model.geo.addLine(48+offset, 45+offset, 48+offset)
        gmsh.model.geo.addLine(41+offset, 45+offset, 49+offset)
        gmsh.model.geo.addLine(44+offset, 48+offset, 52+offset)
    else:
        gmsh.model.geo.addLine(42+offset2, 42+offset, 41+offset)
        gmsh.model.geo.addLine(43+offset, 43+offset2, 43+offset)
        gmsh.model.geo.addLine(46+offset2, 46+offset, 45+offset)
        gmsh.model.geo.addLine(47+offset, 47+offset2, 47+offset)
        gmsh.model.geo.addLine(21+offsetl, 46+offset2, 53+offset)
        gmsh.model.geo.addLine(24+offsetl, 47+offset2, 56+offset)
        gmsh.model.geo.addLine(25+offset, 42+offset2, 61+offset)
        gmsh.model.geo.addLine(28+offset, 43+offset2, 64+offset)
    
    gmsh.model.geo.addLine(42+offset, 43+offset, 42+offset)
    gmsh.model.geo.addLine(46+offset, 47+offset, 46+offset)
    gmsh.model.geo.addLine(42+offset, 46+offset, 50+offset)
    gmsh.model.geo.addLine(43+offset, 47+offset, 51+offset)
    
    gmsh.model.geo.addLine(22+offsetl, 46+offset, 54+offset)
    gmsh.model.geo.addLine(23+offsetl, 47+offset, 55+offset)
    gmsh.model.geo.addLine(26+offset, 42+offset, 62+offset)
    gmsh.model.geo.addLine(27+offset, 43+offset, 63+offset)

    gmsh.model.geo.addLine(25+offset, 21+offsetl, 65+offset)
    gmsh.model.geo.addLine(26+offset, 22+offsetl, 66+offset)
    gmsh.model.geo.addLine(27+offset, 23+offsetl, 67+offset)
    gmsh.model.geo.addLine(28+offset, 24+offsetl, 68+offset)
    

    if first:
        gmsh.model.geo.addCurveLoop([41+offset, 50+offset, -45-offset, -49-offset2], 40+offset)
        gmsh.model.geo.addCurveLoop([43+offset, 52+offset2, -47-offset, -51-offset], 42+offset)
        gmsh.model.geo.addCurveLoop([-53-offset, -24-offsetl, 56+offset, 48+offset2], 47+offset)
        gmsh.model.geo.addCurveLoop([-61-offset, -28-offset, 64+offset, 44+offset2], 52+offset)
        gmsh.model.geo.addCurveLoop([65+offset, 53+offset, -49-offset2, -61-offset], 53+offset)
        gmsh.model.geo.addCurveLoop([68+offset, 56+offset, -52-offset2, -64-offset], 56+offset)

        gmsh.model.geo.addCurveLoop([44+offset, 49+offset, -48-offset, -52-offset], 43+offset)
    else:
        gmsh.model.geo.addCurveLoop([41+offset, 50+offset, -45-offset, -50-offset2], 40+offset)
        gmsh.model.geo.addCurveLoop([43+offset, 51+offset2, -47-offset, -51-offset], 42+offset)
        gmsh.model.geo.addCurveLoop([-53-offset, -24-offsetl, 56+offset, -46-offset2], 47+offset)
        gmsh.model.geo.addCurveLoop([-61-offset, -28-offset, 64+offset, -42-offset2], 52+offset)
        gmsh.model.geo.addCurveLoop([65+offset, 53+offset, -50-offset2, -61-offset], 53+offset)
        gmsh.model.geo.addCurveLoop([68+offset, 56+offset, -51-offset2, -64-offset], 56+offset)

    gmsh.model.geo.addCurveLoop([42+offset, 51+offset, -46-offset, -50-offset], 41+offset)

    # gmsh.model.geo.addCurveLoop([41+offset, 42+offset, 43+offset, 44+offset2], 44+offset)
    # gmsh.model.geo.addCurveLoop([45+offset, 46+offset, 47+offset, 48+offset2], 45+offset)

    gmsh.model.geo.addCurveLoop([45+offset, -54-offset, -21-offsetl, 53+offset], 44+offset)
    gmsh.model.geo.addCurveLoop([54+offset, 46+offset, -55-offset, -22-offsetl], 45+offset)
    gmsh.model.geo.addCurveLoop([-23-offsetl, 55+offset, 47+offset, -56-offset], 46+offset)

    gmsh.model.geo.addCurveLoop([41+offset, -62-offset, -25-offset, 61+offset], 49+offset)
    gmsh.model.geo.addCurveLoop([62+offset, 42+offset, -63-offset, -26-offset], 50+offset)
    gmsh.model.geo.addCurveLoop([-27-offset, 63+offset, 43+offset, -64-offset], 51+offset)

    gmsh.model.geo.addCurveLoop([66+offset, 54+offset, -50-offset, -62-offset], 54+offset)
    gmsh.model.geo.addCurveLoop([67+offset, 55+offset, -51-offset, -63-offset], 55+offset)

    gmsh.model.geo.addCurveLoop([65+offset, -24-offsetl, -68-offset, 28+offset], 57+offset)
    gmsh.model.geo.addCurveLoop([66+offset, 22+offsetl, -67-offset, -26-offset], 58+offset)
    gmsh.model.geo.addCurveLoop([67+offset, 23+offsetl, -68-offset, -27-offset], 59+offset)
    gmsh.model.geo.addCurveLoop([66+offset, -21-offsetl, -65-offset, 25+offset], 60+offset)
    
    # Surfaces stiffner
    for su in range(40 + offset, 43 + offset):
        gmsh.model.geo.addPlaneSurface([su], su)
    if first:
        gmsh.model.geo.addPlaneSurface([43 + offset], 43 + offset)
    # gmsh.model.geo.addPlaneSurface([44 + offset] + curveloops1, 44 + offset)
    # gmsh.model.geo.addPlaneSurface([45 + offset] + curveloops2, 45 + offset)
    for su in range(44 + offset, 48 + offset):
        gmsh.model.geo.addPlaneSurface([su], su)
    for su in range(49 + offset, 61 + offset):
        gmsh.model.geo.addPlaneSurface([su], su)

    # volumes Stiffner
    gmsh.model.geo.addSurfaceLoop([20 + offsetl, 37 + offset, 57 + offset, 58 + offset, 59 + offset, 60 + offset], 41+offset)
    gmsh.model.geo.addSurfaceLoop([40 + offset, 44 + offset, 49 + offset, 53 + offset, 54 + offset, 60 + offset], 42+offset)
    gmsh.model.geo.addSurfaceLoop([41 + offset, 45 + offset, 50 + offset, 54 + offset, 55 + offset, 58 + offset], 43+offset)
    gmsh.model.geo.addSurfaceLoop([42 + offset, 46 + offset, 51 + offset, 55 + offset, 56 + offset, 59 + offset], 44+offset)
    if first:
        gmsh.model.geo.addSurfaceLoop([43 + offset2, 47 + offset, 52 + offset, 53 + offset, 56 + offset, 57 + offset], 45+offset)
    else:
        gmsh.model.geo.addSurfaceLoop([41 + offset2, 47 + offset, 52 + offset, 53 + offset, 56 + offset, 57 + offset], 45+offset)

    for vo in range(41+offset, 46+offset): 
        gmsh.model.geo.addVolume([vo],vo)

    gmsh.model.geo.synchronize()


# === Main ===

gmsh.initialize()
# gmsh.option.setNumber('General.Verbosity', 99)
gmsh.model.add("Assemblage")

meshSizeC = 1e-3
meshSizeG = 1e-3
meshSizeS = 1e-3

CouronnesNbr = 1
PlatesNbr = 3
AnglePlate = 360./PlatesNbr

StiffenerHalfThickness = 2.e-3
TabFuelThickness = [0.6e-3, 0.6e-3, 0.6e-3, 0.6e-3, 0.6e-3, 0.6e-3, 0.6e-3, 0.6e-3]
TabCladThickness = [0.4e-3, 0.4e-3, 0.4e-3, 0.4e-3, 0.4e-3, 0.4e-3, 0.4e-3, 0.4e-3]
TabRayonCintrage = [22.e-3, 25.5e-3, 29.e-3, 32.5e-3, 36.e-3, 39.5e-3, 43.e-3, 46.5e-3]
TabStiffenerIntWidth = [1.5e-3, 1.05e-3, 1.05e-3, 1.05e-3, 1.05e-3, 1.05e-3, 1.05e-3, 1.05e-3]
TabStiffenerExtWidth = [1.05e-3, 1.05e-3, 1.05e-3, 1.05e-3, 1.05e-3, 1.05e-3, 1.05e-3, 1.5e-3]

FuelHeight = 600.e-3
PlateHeight = 700.e-3
StiffenerBottomHeight = 70.e-3
StiffenerTopHeight = 70.e-3

hautGaineB = StiffenerBottomHeight
hautGaineH = hautGaineB + PlateHeight

offsetBase = 100

# PhysicalGroup surfaces

sgh=[]
sge=[]
sgb=[]
sgi=[]
sse=[]
ssb=[]
ssh=[]

# PhysicalGroup volumes

comb=[]
gaine=[]
stiffeners=[]

# Transfinite lines

fuelHeight=[]
fuelLength=[]
fuelThick=[]
cladConn=[]
stifThick=[]
stifConn=[]


# plaques

for i in range(CouronnesNbr):
    for j in range(PlatesNbr):
        offset = (i * PlatesNbr + j) * offsetBase
        FuelThickness = TabFuelThickness[i]
        CladThickness = TabCladThickness[i]
        RayonCintrage = TabRayonCintrage[i]
        angle = j * AnglePlate
        plaque(offset,FuelThickness,CladThickness,RayonCintrage,angle)

        # PhysicalGroup surfaces
        sgb.extend([21+offset])
        sge.extend([22+offset])
        sgh.extend([23+offset])
        sgi.extend([24+offset])

        # PhysicalGroup volumes
        comb.extend([1+offset])
        gaine.extend([21+offset, 22+offset, 23+offset, 24+offset, 25+offset, 26+offset])

        # groupe pour transfinite
        fuelHeight.extend([2+offset, 4+offset, 6+offset, 8+offset])
        fuelHeight.extend([22+offset, 24+offset, 26+offset, 28+offset])
        fuelLength.extend([9+offset, 10+offset, 11+offset, 12+offset])
        fuelLength.extend([29+offset, 30+offset, 31+offset, 32+offset])
        fuelThick.extend([1+offset, 3+offset, 5+offset, 7+offset])
        fuelThick.extend([21+offset, 23+offset, 25+offset, 27+offset])
        cladConn.extend([33+offset, 34+offset, 35+offset, 36+offset])
        cladConn.extend([37+offset, 38+offset, 39+offset, 40+offset])


# stiffeners 

hautStifB = 0.
hautStifH = hautGaineH + StiffenerTopHeight

for i in range(CouronnesNbr):
    for j in range(PlatesNbr):
        rintStif = TabRayonCintrage[i] - TabFuelThickness[i] / 2. - TabCladThickness[i] - TabStiffenerIntWidth[i]
        rextStif = TabRayonCintrage[i] + TabFuelThickness[i] / 2. + TabCladThickness[i] + TabStiffenerExtWidth[i]

        offset = (i * PlatesNbr + j) * offsetBase
        angStif = 60. + j * AnglePlate
        stiffener(offset, rintStif, rextStif, hautStifB, hautStifH, angStif, i, j)

        # PhysicalGroup surfaces
        sse.extend([44+offset,45+offset,46+offset,47+offset])
        sse.extend([49+offset,50+offset,51+offset,52+offset])
        ssb.extend([40+offset])
        ssh.extend([42+offset])

        # PhysicalGroup volumes
        stiffeners.extend([41+offset,42+offset,43+offset,44+offset,45+offset])
        
        # groupe pour transfinite
        fuelThick.extend([41+offset, 43+offset, 45+offset, 47+offset])
        if (i == 0):
            fuelHeight.extend([42+offset, 44+offset, 46+offset, 48+offset])
            stifThick.extend([49+offset, 50+offset, 51+offset, 52+offset])
        else:
            fuelHeight.extend([42+offset, 46+offset])
            stifThick.extend([50+offset, 51+offset])
        stifThick.extend([65+offset, 66+offset, 67+offset, 68+offset])
        stifConn.extend([53+offset, 54+offset, 55+offset, 56+offset])
        stifConn.extend([61+offset, 62+offset, 63+offset, 64+offset])

# Synchronisation finale
gmsh.model.geo.synchronize()


# PhysicalGroup volumes

gmsh.model.addPhysicalGroup(3, comb, 1) 
gmsh.model.setPhysicalName(3, 1, "comb")

gmsh.model.addPhysicalGroup(3, gaine, 2) 
gmsh.model.setPhysicalName(3, 2, "gaine")

gmsh.model.addPhysicalGroup(3, stiffeners, 3)  
gmsh.model.setPhysicalName(3, 3, "stiffeners")


# PhysicalGroup surfaces

gmsh.model.addPhysicalGroup(2, sgh, 4)
gmsh.model.setPhysicalName(2, 4, "sgh")

gmsh.model.addPhysicalGroup(2, sge, 5)
gmsh.model.setPhysicalName(2, 5, "sge")

gmsh.model.addPhysicalGroup(2, sgb, 6)
gmsh.model.setPhysicalName(2, 6, "sgb")

gmsh.model.addPhysicalGroup(2, sgi, 7)
gmsh.model.setPhysicalName(2, 7, "sgi")

gmsh.model.addPhysicalGroup(2, sse, 8) 
gmsh.model.setPhysicalName(2, 8, "sse")

gmsh.model.addPhysicalGroup(2, ssb, 9)
gmsh.model.setPhysicalName(2, 9, "ssb")

gmsh.model.addPhysicalGroup(2, ssh, 10)
gmsh.model.setPhysicalName(2, 10, "ssh")


# Maillage transfinis et recombinaison automatique
#gmsh.model.mesh.setTransfiniteAutomatic(cornerAngle=10, recombine=True)  # 3.06177 rad = 171,887 degré

# Création du parser d'arguments
parser = argparse.ArgumentParser(description="Script de génération de maillage")

parser.add_argument("--densHaut",       type=int, default=10, help="Nombre d'éléments dans la hauteur")
parser.add_argument("--densFuelLength", type=int, default=15, help="Nombre d'éléments dans le long du combustible")
parser.add_argument("--densFuelThick",  type=int, default=5, help="Nombre d'éléments dans l'épaisseur du combustible")
parser.add_argument("--densCladConn",   type=int, default=5, help="Nombre d'éléments dans l'épaisseur de la gaine")
parser.add_argument("--densStifThick",  type=int, default=5, help="Nombre d'éléments dans l'épaisseur du raidisseur")
parser.add_argument("--densStifConn",   type=int, default=5, help="Nombre d'éléments dans le raidisseur hors gaine")
parser.add_argument("--output_file",    type=str, default="assemblage_hexa.msh", help="Chemin du fichier de sortie")

# Lecture des arguments
args = parser.parse_args()

# Utilisation des valeurs passées en ligne de commande ou valeurs par défaut
densHaut = args.densHaut
densFuelLength = args.densFuelLength
densFuelThick = args.densFuelThick
densCladConn = args.densCladConn
densStifThick = args.densStifThick
densStifConn = args.densStifConn

for line in fuelHeight:
    gmsh.model.mesh.setTransfiniteCurve(line,densHaut)
for line in fuelLength:
    gmsh.model.mesh.setTransfiniteCurve(line,densFuelLength)
for line in fuelThick:
    gmsh.model.mesh.setTransfiniteCurve(line,densFuelThick)
for line in cladConn:
    gmsh.model.mesh.setTransfiniteCurve(line,densCladConn)
for line in stifThick:
    gmsh.model.mesh.setTransfiniteCurve(line,densStifThick)
for line in stifConn:
    gmsh.model.mesh.setTransfiniteCurve(line,densStifConn)

for s in gmsh.model.getEntities(2):
    gmsh.model.mesh.setTransfiniteSurface(s[1])
    gmsh.model.mesh.setRecombine(s[0], s[1])
for v in gmsh.model.getEntities(3):
    gmsh.model.mesh.setTransfiniteVolume(v[1])


# Options pour le maillage

# gmsh.option.setNumber("Mesh.Optimize", 1)  # Optimiser le maillage
# gmsh.option.setNumber("Mesh.Smoothing", 1)  # Appliquer un lissage au maillage

gmsh.option.setNumber("Mesh.ElementOrder", 1)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.Format", 0)  

t0 = time.time()

# Génération du maillage 
gmsh.model.mesh.generate(1)  # mesh1 1D
t1 = time.time()
print("Timing 1D=", t1 - t0, "seconds")

gmsh.model.mesh.generate(2)  # mesh2 2D
t2 = time.time()
print("Timing 2D=", t2 - t1, "seconds")

gmsh.model.mesh.generate(3)  # mesh3 3D
t3 = time.time()
print("Timing 3D=", t3 - t2, "seconds")

#Sauvegarder le maillage dans un fichier
gmsh.write(args.output_file)
print(f"Fichier généré : {args.output_file}")
print(f"Temps génération maillage : {time.time() - t0:.2f} s")

# Finaliser
gmsh.finalize()
