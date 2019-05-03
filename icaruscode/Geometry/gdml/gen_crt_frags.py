'''Generate the GDML for the SBND CRT geometry.

We define a global XML DOM which the function for each volume will append to.
Volume-building functions are called hierarchically, so if you call module(),
it will construct the module and all the parts that make it up, so you end
up with complete GDML for one module.

Each physical volume has a corresponding unique logical volume, as required
by LArG4 to keep track of energy depositions. The solids, however, can safely
be referenced many times, and so are stored only once (using a hash keyed on
the the linear dimensions).

The output of this code is a file "crt.gdml" which contains the GDML snippets
to paste into the full SBND geometry.

Created by A. Mastbaum <mastbaum@uchicago.edu>, 2016/10/27
Downloaded by C. Hilgenberg <chilge@rams.colostate.edu> 2017/10/25
Modified by C. Hilgenberg 2017/10/26
'''

import math
import xml.etree.cElementTree as ET
from xml.dom import minidom

#################### Parameters #####################
#warm vessel (cm)
WVWIDTH  = 972.0 
WVHEIGHT = 614.0
WVLENGTH = 2209.0
ISLANDWIDTH = 118.0

#strip width, length, thickness
YM = 4.1   #MINOS
ZM = 800.0
XM = 1.0

XC = 23.0   #CERN
ZC = 184.0
YC = 1.5

XD = 5.0    #Double Chooz
ZD = 322.5
YD = 1.0

#Strips per layer
NXM = 20
NXC = 8
NXD = 32

# cm padding between strips and module (Al thickness)
PADM = 0.05  
PADC = 0.1
PADD = 0.05
PADModule = 0.1
PADStrip = 0.01
PADTagger = 0.001

#MINOS mounting
LAYERSPACE = 10.0 #MINOS edge-to-edge distance between adjacent layers (cm)
NMODSTACK = 9 #number of lateral MINOS modules in a single layer, single stack
SLIDERSPACE = 25.0 #MINOS edge-to-edge distance between fixed and sliding stacks (cm)
STACKOVERLAP = 50.0 #MINOS stack horizontal overlap (cm)

#DC mounting
DCSPACER=32.6 #foam spacer between DC modules in rows of 5 (strip normal to drift direction) (cm)
LONGOFF5=(3*ISLANDWIDTH+481.8)*0.5+181.8
LONGOFF2= (ISLANDWIDTH+181.8)*0.5

#CERN mounting
NTOPX=6
NTOPZ=14
NSLOPELAT=13
NSLOPEFRONT=6
SLOPEINCLINATION=60.0 #degrees w.r.t. vertical

#MINOS sections positions
posMINOSSide1InDetEncl = (-537.13, -81.735, 0)
posMINOSSide2InDetEncl = (537.13, -81.735, 0)
posMINOSFrontInDetEncl = (0,-51.1349999999999, -1173.6)
posMINOSBackInDetEncl= (0,-51.1349999999999, 1173.6)

#DC section positions
posDCInDetEncl = (0,-480.135, 0)

#CERN sections positions
posCERNTopInDetEncl = (0, 479.565, 0)
posCERNFrontInDetEncl = (0, 397.0, -1306.84)
posCERNBackInDetEncl = (0, 397.0, 1306.84)
posCERNLeftInDetEncl = (570.04, 397.0, 0)
posCERNRightInDetEncl = (-570.04, 397.0, 0)


########################################################

gdml = ET.Element('gdml')
#materials = ET.SubElement(gdml, 'materials')
solids = ET.SubElement(gdml, 'solids')
structure = ET.SubElement(gdml, 'structure')
solids_store = {}

mod_id = -1
nModM = 0
nModC = 0
nModD = 0

def get_mod_id(style='m'):
    global mod_id
    global nModM
    global nModC
    global nModD

    mod_id += 1
    if style == 'm':
        nModM += 1
    if style == 'c':
        nModC += 1
    if style == 'd':
        nModD += 1
    return mod_id

def get_mod_id_num():
    global mod_id
    return str(mod_id)

def strip(style="m", modnum=0, stripnum=0):
    '''Build one scintillator strip.'''

    if style=="m":
        x=XM
        y=YM
        z=ZM
        name = 'MINOS'
    if style=="c":
        x=XC
        y=YC
        z=ZC
        name = 'CERN'
    if style=="d":
        x=XD
        y=YD
        z=ZD
        name = 'DC'

    xx = str(x)
    yy = str(y)
    zz = str(z) 

    sname = 'AuxDetSensitive_' + name + '_strip'
    vname = 'volAuxDetSensitive_'
    vname += name
    vname += '_module_'

    if modnum < 10:
        vname += '00'
    elif modnum < 100: 
        vname += '0'

    vname += str(modnum) + '_strip_'

    if stripnum < 10: 
        vname += '0'
    vname += str(stripnum)

    if not sname in solids_store:
        s = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz) #g4 solid
        solids_store[sname] = s

    else:
        s = solids_store[sname]

    #vname = 'volAuxDetSensitive_' + name
    v = ET.SubElement(structure, 'volume', name=vname) #Logical volume
    ET.SubElement(v, 'materialref', ref='Polystyrene')
    ET.SubElement(v, 'solidref', ref=sname)

    #print("strip produced!")

    return s, v #return solid, logical volumes

def module(style="m", reg='tt'):
    '''Build an edge-to-edge array of scintillator strips.'''

    if style=="m":
        x=XM
        y=YM
        z=ZM
        ny=NXM
        xx = str(x+2*PADM+2*PADStrip)
        yy = str(y*ny+2*PADM+(ny+1)*PADStrip)
        zz = str(z+2*PADM+2*PADStrip)
        xxsub = str(x+2*PADStrip)
        yysub = str(y*ny+(ny+1)*PADStrip)
        zzsub = str(z+2*PADStrip)
	name = "MINOS"

    if style=="c":
        x=XC
        y=YC
        z=ZC
        ny=NXC
        xx = str(x*ny+2*PADC+(ny+1)*PADStrip)
        yy = str(y*2+2*PADC+3*PADStrip)
        zz = xx
        xxsub = str(x*ny+(ny+1)*PADStrip)
        yysub = str(y*2+3*PADStrip)
        zzsub = xxsub
	name = "CERN"

    if style=="d":
        x=XD
        y=YD
        z=ZD
        ny=NXD
        xx = str(x*(ny+0.5)+2*PADD+(ny+2)*PADStrip)
        yy = str(y*2+2*PADD+3*PADStrip)
        zz = str(z+2*PADD+2*PADStrip)
        xxsub = str(x*(ny+0.5)+(ny+2)*PADStrip)
        yysub = str(y*2+3*PADStrip)
        zzsub = str(z+2*PADStrip)
	name = "DC"

    modnum = get_mod_id(style)
    stripnum = 0

    sname = 'AuxDet_' + name + '_module'
    vname = 'vol' + sname + '_'

    if modnum < 10:
        vname += '00'
    elif modnum < 100:
        vname += '0'
    vname += str(modnum)
    vname +='_'

    if reg=='tt': vname += 'Top'
    if reg=='sf': vname += 'SlopeFront'        
    if reg=='sb': vname += 'SlopeBack'
    if reg=='sl': vname += 'SlopeLeft'
    if reg=='sr': vname += 'SlopeRight'        
    if reg=='ff': vname += 'Front'
    if reg=='bb': vname += 'Back'
    if reg=='ll': vname += 'Left'        
    if reg=='rr': vname += 'Right'
    if reg=='bt': vname += 'Bottom'

    snamein  = sname+'_inner'

    if not sname in solids_store: 

        s = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
        sin = ET.SubElement(solids, 'box', name=snamein, lunit="cm", x=xxsub, y=yysub, z=zzsub)

        solids_store[sname] = s
        solids_store[snamein] = sin

    else:
        s = solids_store[sname]
        sin = solids_store[snamein]

    strips  = []
    strips2 = []

    for i in range(ny):
        strips.append(strip(style, modnum, stripnum))
        stripnum += 1

    if style=='d' or style=='c':    
        for i in range(ny):
            strips2.append(strip(style, modnum, stripnum))
            stripnum += 1

    vnamein = vname + '_inner'
    vin = ET.SubElement(structure, 'volume', name=vnamein)
    ET.SubElement(vin, 'materialref', ref='Air')
    ET.SubElement(vin, 'solidref', ref=snamein)

    for i, (es, ev) in enumerate(strips):
        pv = ET.SubElement(vin, 'physvol')
        ET.SubElement(pv, 'volumeref', ref=ev.attrib['name'])

        if style=='m':
            dy = (2*i - ny + 1)* 0.5 * (y+PADStrip)
            dx=0
        if style=='c':
            dx = (2*i - ny + 1)* 0.5 * (x+PADStrip)
            dy=0.5*(y+PADStrip)

        if style=='d':
            dy= 0.5*(y+PADStrip)
            dx=(i - 0.5*ny + 0.25) * (x+PADStrip)

        posname = 'pos' + ev.attrib['name']
        ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(dx), y=str(dy), z='0')

    if style=='c':
        for i, (es, ev) in enumerate(strips2):
            pv = ET.SubElement(vin, 'physvol')
            ET.SubElement(pv, 'volumeref', ref=ev.attrib['name'])

            dy= -0.5*(y+PADStrip)
            dz=(2*i - ny + 1)* 0.5 * (x+PADStrip)

            posname = 'pos' + ev.attrib['name']
            ET.SubElement(pv, 'position', name=posname,
                          unit="cm", x='0', y=str(dy), z=str(dz))
            posname = 'rot' + ev.attrib['name']
            ET.SubElement(pv, 'rotation', name=posname, unit="deg", x='0', y='90', z='0')


    if style=='d':
        for i, (es, ev) in enumerate(strips2):
            pv = ET.SubElement(vin, 'physvol')
            ET.SubElement(pv, 'volumeref', ref=ev.attrib['name'])

            dy= -0.5*(y+PADStrip)
            dx=(i - 0.5*ny + 0.75) * (x+PADStrip)

            posname = 'pos' + ev.attrib['name']
            ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(dx), y=str(dy), z='0')
            #DC strips centered (FIX ME!!)

    v = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(v, 'materialref', ref='ALUMINUM_Al')
    ET.SubElement(v, 'solidref', ref=sname)

    pv = ET.SubElement(v, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vin.attrib['name'])

    #print ('module produced!')

    return s, v

def minosSideTagger(side='L', x0=0, y0=0, z0=0):
    ''' Build a side tagger (3 stacks)
    '''
    coords = []
    modules = []

    xx = str(SLIDERSPACE+2*LAYERSPACE+4*(XM+2*(PADM+PADStrip))+2*PADTagger)
    yy = str(NMODSTACK * (YM*NXM+2*PADM+(NXM+1)*PADStrip) + (NMODSTACK+1)*PADModule + 2*PADTagger)
    zz = str(3 * (ZM + 2*(PADM+PADStrip)) - 2*STACKOVERLAP + 2*PADTagger)

	#loop over stacks
    for layer in range (6):
        if (layer==0 or layer==1):
            dz=-1*(ZM + 2*(PADM+PADStrip) - STACKOVERLAP)
            dx=-1*((SLIDERSPACE+LAYERSPACE)*0.5 + XM + 2*(PADM+PADStrip))
        if (layer==2 or layer==3):
            dz=0
            dx=(SLIDERSPACE+LAYERSPACE)*0.5 + XM + 2*(PADM+PADStrip)
        if (layer==4 or layer==5):
            dz=ZM + 2*(PADM+PADStrip) - STACKOVERLAP
            dx=-1*((SLIDERSPACE+LAYERSPACE)*0.5 + XM + 2*(PADM+PADStrip))
        if (side=='L'): dx *= -1

        dx+=((-1)**layer)*(LAYERSPACE+XM+2*(PADM+PADStrip))/2.0

        #loop over modules in stack
        for i in range(NMODSTACK):

            dy = (YM*NXM + 2*PADM+(NXM+1)*PADStrip+PADModule) * (0.5 * (1-NMODSTACK) + i )

            coords.append((x0+dx,y0+dy,z0+dz))

    for i in range(len(coords)):
        if side=='L':
            modules.append(module('m','ll'))
        if side=='R':
            modules.append(module('m','rr'))

    sname = 'tagger_'
    if side == 'L':
        sname += 'left'
    if side == 'R':
        sname += 'Right'
    vname = 'vol_'+ sname

    stagger = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
    vtagger = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(vtagger, 'materialref', ref='Air')
    ET.SubElement(vtagger, 'solidref', ref=sname)

    #place left side module phy. vol.s
    for i, (xc,yc,zc) in enumerate(coords):

        (s,v)=modules[i]
        pv = ET.SubElement(vtagger, 'physvol')
        ET.SubElement(pv, 'volumeref', ref=v.attrib['name'])

        posname = 'pos' + v.attrib['name']
        ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    return stagger, vtagger

def minosFrontTagger(side='U', x0=0, y0=0, z0=0):
    ''' Build front MINOS tagger (2 layers in X-Y) on upstream face
    '''
    nmody = 11
    coords = []
    modules = []

    xx = str( max ( nmody*(YM*NXM+2*PADM+(NXM+1)*PADStrip)+(nmody+1)*PADModule, ZM+2*(PADM+PADStrip)) + 2*PADTagger)
    yy = str( max ( NMODSTACK*(YM*NXM+2*PADM+(NXM+1)*PADStrip)+(NMODSTACK-1)*PADModule, ZM+2*(PADM+PADStrip))+2*PADTagger)
    zz = str( 2*(XM+2*(PADM+PADStrip)) + LAYERSPACE + 2*PADTagger )

    dz = (LAYERSPACE + XM + 2*(PADM+PADStrip))*0.5
    if side=='D' : dz*=-1

    for i in range(nmody):

        dx = (YM*NXM + 2*PADM + (NXM+1)*PADStrip + PADModule ) * (0.5 * (1-11) + i )

        coords.append((dx+x0,y0,dz+z0,1)) #x,y,z,vert=true

    dz*=-1

    for i in range(NMODSTACK):

        dy = (YM*NXM + 2*PADM + (NXM+1)*PADStrip + PADModule ) * (0.5 * (1-NMODSTACK) + i )
    
        coords.append((x0,dy+y0,dz+z0,0)) #x,y,z,vert=false

    for i in range(len(coords)):
        if side=='U':
            modules.append(module('m','ff'))
        if side=='D':
            modules.append(module('m','bb'))

    sname = 'tagger_'
    if side == 'U':
        sname += 'Front'
    if side == 'D':
        sname += 'Back'
    vname = 'vol_'+ sname

    stagger = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
    vtagger = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(vtagger, 'materialref', ref='Air')
    ET.SubElement(vtagger, 'solidref', ref=sname)

    #place left side module phy. vol.s
    for i, (xc,yc,zc,r) in enumerate(coords):

        (s,v)=modules[i]
        pv = ET.SubElement(vtagger, 'physvol')
        ET.SubElement(pv, 'volumeref', ref=v.attrib['name'])

        posname = 'pos' + v.attrib['name']
        ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

        posname = 'rot' + v.attrib['name']
        if r==1 : ET.SubElement(pv, 'rotation', name=posname, unit="deg", x='90', y='0', z='90')
        if r==0 : ET.SubElement(pv, 'rotation', name=posname, unit="deg", x='0', y='90', z='0')

    return stagger, vtagger


def DCTagger(x0=0, y0=0, z0=0):
    ''' Build bottom tagger
    '''
    modwidth = XD*(NXD+0.5)+2*PADD+(NXD+2)*PADStrip
    xx = str(modwidth*5 + DCSPACER*4 + 2*PADTagger)
    yy = str(2*(YD+PADD+PADTagger)+3*PADStrip)
    zz = str(WVLENGTH)

    coords = []
    modules = []
    rot = 0
    for i in range(14):

        if (i<5):
            dx = (2*i-5+1)*0.5*(modwidth+ DCSPACER)
            dz = -1*LONGOFF5
        if (i==5 or i==6):
            dx = (ZD + 2*(PADD+PADStrip))*0.5*(-1)**i
            dz = -1*LONGOFF2
        if (i==7 or i==8):
            dx = (ZD + 2*(PADD+PADStrip))*0.5*(-1)**i
            dz = LONGOFF2
        if (i>8):
            dx = (2*(i-9)-5+1)*0.5*(modwidth+ DCSPACER)
            dz = LONGOFF5

        if (i>4 and i<9):
            rot = 1
        else :
            rot = 0

        coords.append((dx+x0,y0,dz+z0,rot))

    for i in range(len(coords)):
        modules.append(module('d','bt'))

    sname = 'tagger_Bottom'
    vname = 'vol_'+ sname

    stagger = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
    vtagger = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(vtagger, 'materialref', ref='Air')
    ET.SubElement(vtagger, 'solidref', ref=sname)

    #place left side module phy. vol.s
    for i, (xc,yc,zc,r) in enumerate(coords):

        (s,v)=modules[i]
        pv = ET.SubElement(vtagger, 'physvol')
        ET.SubElement(pv, 'volumeref', ref=v.attrib['name'])

        posname = 'pos' + v.attrib['name']
        ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

        if r==1 :
            posname = 'rot' + v.attrib['name']
            ET.SubElement(pv, 'rotation', name=posname, unit="deg", x='0', y='90', z='0')

    return stagger, vtagger


def cernTopTagger(x0=0, y0=0, z0=0):
    ''' Build front MINOS tagger (2 layers in X-Y) on upstream face
    '''
    modwidth = ZC + 2*PADC + (NXC+1)*PADStrip
    xx = str(NTOPX*modwidth+2*PADTagger+(NTOPX-1)*PADModule)
    yy = str(2*YC+3*PADStrip+2*PADC+2*PADTagger)
    zz = str(NTOPZ*modwidth + 2*PADTagger + (NTOPZ-1)*PADModule)

    coords = []
    modules = []

    dz = 0.5*(modwidth+PADModule)*(1 - NTOPZ)
    dx = 0.5*(modwidth+PADModule)*(1 - NTOPX)

    for i in range(NTOPX*NTOPZ):

        coords.append((dx+x0,y0,dz+z0))

        if (i+1)%NTOPZ == 0:
            dx+= modwidth + PADModule
            dz = 0.5*(modwidth+PADModule)*(1 - NTOPZ)
        else: dz+= modwidth + PADModule
    
    for i in range(len(coords)):
        modules.append(module('c','tt'))

    sname = 'tagger_Top'
    vname = 'vol_'+ sname

    stagger = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
    vtagger = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(vtagger, 'materialref', ref='Air')
    ET.SubElement(vtagger, 'solidref', ref=sname)

    #place left side module phy. vol.s
    for i, (xc,yc,zc) in enumerate(coords):

        (s,v)=modules[i]
        pv = ET.SubElement(vtagger, 'physvol')
        ET.SubElement(pv, 'volumeref', ref=v.attrib['name'])

        posname = 'pos' + v.attrib['name']
        ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    return stagger, vtagger

def cernSlopeSideTagger(side='L',x0=0, y0=0, z0=0):
    ''' Build sloped CERN tagger
    '''
    modwidth = ZC + 2*PADC + (NXC+1)*PADStrip
    xx = str(modwidth+2*PADTagger)
    yy = str(2*YC+3*PADStrip+2*PADC+2*PADTagger)
    zz = str(NSLOPELAT*modwidth + 2*PADTagger + (NSLOPELAT-1)*PADModule)

    coords = []
    modules = []

    dz = 0.5*(modwidth+PADModule)*(1 - NSLOPELAT)

    for i in range(NSLOPELAT):

        coords.append((x0,y0,dz+z0))

        dz+= modwidth+PADModule
    
    for i in range(len(coords)):
        if side == 'L':
            modules.append(module('c','sl'))
        if side == 'R':
            modules.append(module('c','sr'))

    sname = 'tagger_'
    if side == 'L':
        sname += 'SlopeLeft'
    if side == 'R':
        sname += 'SlopeRight'
    vname = 'vol_'+ sname

    stagger = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
    vtagger = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(vtagger, 'materialref', ref='Air')
    ET.SubElement(vtagger, 'solidref', ref=sname)

    #place left side module phy. vol.s
    for i, (xc,yc,zc) in enumerate(coords):

        (s,v)=modules[i]
        pv = ET.SubElement(vtagger, 'physvol')
        ET.SubElement(pv, 'volumeref', ref=v.attrib['name'])

        posname = 'pos' + v.attrib['name']
        ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    return stagger, vtagger


def cernSlopeFrontTagger(side='U',x0=0, y0=0, z0=0):
    ''' Build sloped CERN tagger
    '''
    modwidth = ZC + 2*PADC + (NXC+1)*PADStrip
    xx = str(NSLOPEFRONT*modwidth+2*PADTagger+(NSLOPEFRONT-1)*PADModule)
    yy = str(2*YC+3*PADStrip+2*PADC+2*PADTagger)
    zz = str(modwidth + 2*PADTagger)

    coords = []
    modules = []

    dx = 0.5*(modwidth+PADModule)*(1 - NSLOPEFRONT)

    for i in range(NSLOPEFRONT):

        coords.append((dx+x0,y0,z0))

        dx+= modwidth+PADModule
    
    for i in range(len(coords)):
        if side == 'U':
            modules.append(module('c','sf'))
        if side == 'D':
            modules.append(module('c','sb'))

    sname = 'tagger_'
    if side == 'U':
        sname += 'SlopeFront'
    if side == 'D':
        sname += 'SlopeBack'
    vname = 'vol_'+ sname

    stagger = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
    vtagger = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(vtagger, 'materialref', ref='Air')
    ET.SubElement(vtagger, 'solidref', ref=sname)

    #place left side module phy. vol.s
    for i, (xc,yc,zc) in enumerate(coords):

        (s,v)=modules[i]
        pv = ET.SubElement(vtagger, 'physvol')
        ET.SubElement(pv, 'volumeref', ref=v.attrib['name'])

        posname = 'pos' + v.attrib['name']
        ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    return stagger, vtagger

def detectorEnclosure():
    xx='1235.96'
    yy='963.37'
    zz='2709.56'

    xxint = '1024.578'
    yyint = '955.698'
    zzint = '2334.758'

    xxext = '1235.225'
    yyext = '962.602'
    zzext = '2708.825'

    (s,vll) = minosSideTagger('L',0,0,0) #MINOS Left
    (s,vrr) = minosSideTagger('R',0,0,0) #MINOS Right
    (s,vff) = minosFrontTagger('U',0,0,0) #MINOS Front
    (s,vbb) = minosFrontTagger('D',0,0,0) #MINOS Back
    (s,vbt) = DCTagger(0,0,0) #DC Bottom
    (s,vtt) = cernTopTagger(0,0,0) #CERN top
    (s,vsl) = cernSlopeSideTagger('L',0,0,0) #CERN SlopeLeft
    (s,vsr) = cernSlopeSideTagger('R',0,0,0) #CERN SlopeRight
    (s,vsf) = cernSlopeFrontTagger('U',0,0,0) #CERN SlopeFront
    (s,vsb) = cernSlopeFrontTagger('D',0,0,0) #CERN SlopeBack

    #DetectorEnclosure
    #sname = 'DetEnclosure'
    sname = 'CRT_Shell'
    snameext = sname+'_external'
    snameint = sname+'_internal'
    sexternal = ET.SubElement(solids, 'box', name=snameext, lunit="cm", x=xxext, y=yyext, z=zzext)
    sinternal = ET.SubElement(solids, 'box', name=snameint, lunit="cm", x=xxint, y=yyint, z=zzint)
    sshell = ET.SubElement(solids, 'subtraction', name=sname)
    ET.SubElement(sshell, 'first', ref=snameext)
    ET.SubElement(sshell, 'second', ref=snameint)

    vname = 'vol'+sname
    vshell = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(vshell, 'materialref', ref='Air')
    ET.SubElement(vshell, 'solidref', ref=sname)

    #position MINOS Left
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vll.attrib['name'])

    xc = posMINOSSide1InDetEncl[0]
    yc = posMINOSSide1InDetEncl[1]
    zc = posMINOSSide1InDetEncl[2]

    posname = 'pos' + vll.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    #Position MINOS Right
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vrr.attrib['name'])

    xc = posMINOSSide2InDetEncl[0]
    yc = posMINOSSide2InDetEncl[1]
    zc = posMINOSSide2InDetEncl[2]

    posname = 'pos' + vrr.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    #position MINOS Front
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vff.attrib['name'])

    xc = posMINOSFrontInDetEncl[0]
    yc = posMINOSFrontInDetEncl[1]
    zc = posMINOSFrontInDetEncl[2]

    posname = 'pos' + vff.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    #position MINOS Back 
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vbb.attrib['name'])

    xc = posMINOSBackInDetEncl[0]
    yc = posMINOSBackInDetEncl[1]
    zc = posMINOSBackInDetEncl[2]

    posname = 'pos' + vbb.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    #position DC Bottom
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vbt.attrib['name'])

    xc = posDCInDetEncl[0]
    yc = posDCInDetEncl[1]
    zc = posDCInDetEncl[2]

    posname = 'pos' + vbt.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))  

    #position CERN Top
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vtt.attrib['name'])

    xc = posCERNTopInDetEncl[0]
    yc = posCERNTopInDetEncl[1]
    zc = posCERNTopInDetEncl[2]

    posname = 'pos' + vtt.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc)) 
   
    #position CERN SlopeLeft
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vsl.attrib['name'])

    xc = posCERNLeftInDetEncl[0]
    yc = posCERNLeftInDetEncl[1]
    zc = posCERNLeftInDetEncl[2]

    posname = 'pos' + vsl.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc)) 

    posname = 'rot' + vsl.attrib['name']
    ET.SubElement(pv, 'rotation', name=posname, unit="deg", x='0', y='0', z=str(SLOPEINCLINATION))
 
    #position CERN SlopeRight
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vsr.attrib['name'])

    xc = posCERNRightInDetEncl[0]
    yc = posCERNRightInDetEncl[1]
    zc = posCERNRightInDetEncl[2]

    posname = 'pos' + vsr.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    posname = 'rot' + vsr.attrib['name']
    ET.SubElement(pv, 'rotation', name=posname, unit="deg", x='0', y='0', z=str(-1*SLOPEINCLINATION))

    #position CERN SlopeFront
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vsf.attrib['name'])

    xc = posCERNFrontInDetEncl[0]
    yc = posCERNFrontInDetEncl[1]
    zc = posCERNFrontInDetEncl[2]

    posname = 'pos' + vsf.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))
    
    posname = 'rot' + vsf.attrib['name']
    ET.SubElement(pv, 'rotation', name=posname, unit="deg", x=str(SLOPEINCLINATION), y='0', z='0')

    #position CERN SlopeBack
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vsb.attrib['name'])

    xc = posCERNBackInDetEncl[0]
    yc = posCERNBackInDetEncl[1]
    zc = posCERNBackInDetEncl[2]

    posname = 'pos' + vsb.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))
   
    posname = 'rot' + vsb.attrib['name']
    ET.SubElement(pv, 'rotation', name=posname, unit="deg", x=str(-1*SLOPEINCLINATION), y='0', z='0')

    return sshell,vshell

#######################################################################

#m = ET.SubElement(materials, 'element', name='aluminum', formula='Al', Z='13')
#ET.SubElement(m, 'atom', value='26.9815')

#m = ET.SubElement(materials, 'element', name='nitrogen', formula='N', Z='7')
#ET.SubElement(m, 'atom', value='14.0067')

#m = ET.SubElement(materials, 'element', name='oxygen', formula='O', Z='8')
#ET.SubElement(m, 'atom', value='15.999')

#m = ET.SubElement(materials, 'element', name='argon', formula='Ar', Z='18')
#ET.SubElement(m, 'atom', value='39.9480')

#m = ET.SubElement(materials, 'element', name='hydrogen', formula='H', Z='1')
#ET.SubElement(m, 'atom', value='1.0079')

#m = ET.SubElement(materials, 'element', name='carbon', formula='C', Z='6')
#ET.SubElement(m, 'atom', value='12.0107')

#m = ET.SubElement(materials, 'material', name='ALUMINUM_Al', formula='ALUMINUM_Al')
#ET.SubElement(m, 'D', value='2.6990', unit='g/cm3')
#ET.SubElement(m, 'fraction', n='1.000', ref='aluminum')

#m = ET.SubElement(materials, 'material', name='Air')
#ET.SubElement(m, 'D', value='0.001205', unit='g/cm3')
#ET.SubElement(m, 'fraction', n='0.781154', ref='nitrogen')
#ET.SubElement(m, 'fraction', n='0.209476', ref='oxygen')
#ET.SubElement(m, 'fraction', n='0.00934', ref='argon') 

#m = ET.SubElement(materials, 'material', name='Polystyrene')
#ET.SubElement(m, 'D', value='1.19', unit='g/cm3')
#ET.SubElement(m, 'fraction', n='0.077418', ref='hydrogen')
#ET.SubElement(m, 'fraction', n='0.922582', ref='carbon')

#ws = ET.SubElement(solids, 'box', name='World', lunit="cm", x='10000', y='10000', z='10000')
#w = ET.SubElement(structure, 'volume', name='volWorld')
#ET.SubElement(w, 'materialref', ref='Air')
#ET.SubElement(w, 'solidref', ref='World')

(s,v) = detectorEnclosure()
#(s,v) = DCTagger(0,0,0)
#(s,v) = module('d','bb')
#(s,v) = strip('m')
#ws = ET.SubElement(solids, 'box', name='World', lunit="cm", x='3000', y='3000', z='3000')
#w = ET.SubElement(structure, 'volume', name='volWorld')
#ET.SubElement(w, 'materialref', ref='Air')
#ET.SubElement(w, 'solidref', ref='World')


#pv = ET.SubElement(w, 'physvol')
#ET.SubElement(pv, 'volumeref', ref=v.attrib['name'])

#posname = 'pos' + v.attrib['name']
#ET.SubElement(pv, 'position', name=posname, unit="cm", x='0', y='0', z='0')

print('MINOS modules generated: '+str(nModM))
print('CERN  modules generated: '+str(nModC))
print('DblCh modules generated: '+str(nModD))

# Generate GDML for the world volume, for testing
#ws = ET.SubElement(solids, 'box', name='World', lunit="cm", x='10000', y='10000', z='10000')
#w = ET.SubElement(structure, 'volume', name='volWorld')
#ET.SubElement(w, 'materialref', ref='Air')
#ET.SubElement(w, 'solidref', ref='World')
#pv = ET.SubElement(w, 'physvol')
#ET.SubElement(pv, 'volumeref', ref=tffv.attrib['name'])
#ET.SubElement(pv, 'position', name='posA', unit="cm", x='0', y='0', z='0')
#setup = ET.SubElement(gdml, 'setup', name='Default', version='1.0')
#ET.SubElement(setup, 'world', ref='volWorld')
#mats = ET.parse('mats.gdml')
#gdml.insert(0, mats.getroot())

with open('icarus_crt.gdml', 'w') as f:
    f.write(minidom.parseString(ET.tostring(gdml)).toprettyxml(indent='\t'))

