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
Downloaded by C. Hilgenberg <Chris.Hilgenberg@colostate.edu> 2017/10/25
Modified by C. Hilgenberg 2017/10/26
'''

import csv
import math
import xml.etree.cElementTree as ET
from xml.dom import minidom

#set true to generate standalone CRT shell, false for normal production
testmode = False
printModIds = False

#################### Parameters #####################
#warm vessel (cm), increased as of 7/19/2019 to 
#reflect full extent of WV profile in detector hall
WVWIDTH  = 1031.8  #previously 972.0  
WVHEIGHT = 627.4   #previously 614.0, current value includes support feet, foot height is 30.0cm 
WVLENGTH = 2268.8  #previously 2209.0
WVFOOTELEVATION = 10.16 #height of bottom of WV foot w.r.t. pit floor
ISLANDWIDTH = 118.0 #width of ~square WV support feet islands

IVLENGTH = 1996.0 #length of interior of cold vessel, the inactive + active LAr volume

#OVEROPENZ = 2758.5 - 30.48
TOPLEDGERISERTOFLOOR = 991.9 #pit floor to top of ledge riser upon which the overburden blocks sit
LEDGERISERHEIGTH = 46.99 #wide flange type, W18 x 71
TOPCRTBEAMTOFLOOR = 970.8
CRTBEAMSPACING = 92.71 #horizontal center to center spacing between top CRT support beams

BOTTOMCRTROLLERHEIGHT = 3.02 #distance between bottom CRT module and pit floor
SIDECRTWVOFFSET = 4.13 #set by fiberglass Unistrut standoffs (given as ~4cm, but it's Unitstrut so I assume it was rounded
SIDECRTPOSTWIDTH = 4.13 #Unistruit vertical support posts, dimension normal to side CRT plane
SIDECRTPOSTSPACING = 4.13 #set by Unistruit bracket shelf
SIDECRTSTACKDEPTH = 3*SIDECRTPOSTWIDTH + 2*SIDECRTPOSTSPACING
SIDECRTSHELFTHICK = 0.56 #to be verified, could also be 0.64 depending on steel type

#dimensions of top CRT support beams, wide flange W10 x 49
#true area is larger than that calculated assuming perfect I shape
#91.54 cm^2 vs 92.9 cm^2 
#bottom of CRT beam offset vertically + 0.635cm from bottom of ledge riser
CRTBEAMHEIGHT = 25.3
CRTBEAMWIDTH  = 25.4
CRTBEAMFLANGETHICK = 14.22
CRTBEAMWEBTHICK = 0.86
CRTBEAMMASSDENS = 73 #kg/m, density 7849 kg/m^3

#strip width, length, thickness (cm)
YM = 4.1   #MINOS
ZM = 800.0
XM = 1.0
MINOSSTRAIGHTSNOUT = 37.5
MINOSBENDSNOUT = 26.5
MINOSSNOUTLENGTH = 0.5*(MINOSSTRAIGHTSNOUT+MINOSBENDSNOUT) #guess, to be verified of length along module that snout, where fiber routing occurs, no scintillator, extends


XC = 23.0   #CERN
ZC = 184.0
YCTOP = 1.0
YCBOT = 1.5

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

mModL = ZM+2*PADM+2*PADStrip
mModW = YM*NXM+2*PADM+(NXM+1)*PADStrip
mModH = XM+2*PADM+2*PADStrip
cModW = XC*NXC+2*PADC+(NXC+1)*PADStrip #same as length
cModH = YCTOP+YCBOT+2*PADC+3*PADStrip
dModH = YD*2+2*PADD+3*PADStrip

#MINOS mounting
LAYERSPACE = 8.27 #MINOS edge-to-edge distance between adjacent layers (cm), prevously 10cm
NMODSTACK = 9 #number of lateral MINOS modules in a single layer, single stack
SLIDERSPACE = 18.5+mModH #MINOS center-to-center distance between fixed and sliding stacks' nearest modules (cm), previously 25.0cm
STACKOVERLAP = 50.0 #MINOS stack horizontal overlap (cm)
SIDECRTROLLOFFSET = 44.29 #offset from outmost extend of WV to center fo rolling stack (E/W sides)
SIDECRTNORTHWALLTOFLOOR = 152.2
SIDECRTSOUTHWALLLATOFFSET = 1.1*MINOSSNOUTLENGTH

#DC mounting
DCSPACER=32.6 #foam spacer between DC modules in rows of 5 (strip normal to drift direction) (cm)
LONGOFF5=(3*ISLANDWIDTH+481.8)*0.5+181.8
LONGOFF2= (ISLANDWIDTH+181.8)*0.5

#CERN mounting
CERNMODSPACE = 0.2
NTOPX=6
NTOPZ=14
NSLOPELAT=14
NSLOPEFRONT=6
SLOPEINCLINATION=90.0 #degrees w.r.t. vertical, previously 60 deg
CERNROOFL = NTOPZ*cModW+(NTOPZ-1)*CERNMODSPACE

#CRT shell
SHELLY = 1.1*cModH+TOPCRTBEAMTOFLOOR-BOTTOMCRTROLLERHEIGHT*0.9
#MINOS sections positions
MINOSSOUTHY = -0.5*SHELLY+0.5*(NMODSTACK*mModW+(NMODSTACK-1)*SIDECRTSHELFTHICK+2*PADTagger)+WVFOOTELEVATION+5
MINOSLATFIXY = MINOSSOUTHY
MINOSLATROLLY = MINOSLATFIXY-0.5*mModW+10
MINOSLATSOUTHACTIVEOVERHANG = 2*MINOSSNOUTLENGTH
MINOSLATSOUTHZ = -0.5*WVLENGTH + 0.5*mModL - 0.5*MINOSLATSOUTHACTIVEOVERHANG
MINOSLATNORTHZ = 0.5*WVLENGTH - 0.5*mModL
MINOSLATCENTZ = 0.5*(-0.5*MINOSLATSOUTHACTIVEOVERHANG+MINOSLATSOUTHZ+MINOSLATNORTHZ)

#DC section positions
posDCInDetEncl = (0,-480.135, 0)

#CERN sections positions
CERNRIMSWVOFFSET = 34.0
CERNTOPY = SHELLY*0.5 - 0.6*cModH
CERNRIMSY = CERNTOPY - 0.5*cModH - CRTBEAMHEIGHT - 2.54 - 0.5*cModW
CERNRIMSZ = -0.5*WVLENGTH - CERNRIMSWVOFFSET
CERNRIMNY = CERNRIMSY + 29.0
CERNRIMNZ = CERNROOFL + CERNRIMSZ + 0.6*cModH
CERNTOPZ = CERNRIMSZ + 0.5*CERNROOFL  #roof center assuming south edge of roof aligned with south rim center in z
CERNRIMLATX = 0.5*WVWIDTH + 0.5*cModH + 38.0
CERNRIMLATY = CERNRIMSY
CERNRIMLATZ = CERNTOPZ

SHELLZ = 1.01*(0.5*cModH+CERNRIMNZ-min(CERNRIMSZ+0.5*cModH,MINOSLATSOUTHZ-(mModL+MINOSLATSOUTHACTIVEOVERHANG+SIDECRTSOUTHWALLLATOFFSET)*0.5)) #2*(CERNROOFL - 0.5*IVLENGTH - cModW)
SHELLWVOFFSET = SHELLZ*0.5/1.01 - WVLENGTH*0.5
if CERNRIMSZ+0.5*cModH < MINOSLATSOUTHZ-(mModL+MINOSLATSOUTHACTIVEOVERHANG+SIDECRTSOUTHWALLLATOFFSET)*0.5:
    SHELLWVOFFSET-= CERNRIMSWVOFFSET
else :
    SHELLWVOFFSET-= MINOSLATSOUTHACTIVEOVERHANG
print('SHELL-WV OFFSET: '+str(SHELLWVOFFSET))

#cut MINOS module lengths including snout, index is row number starting from the bottom
minosCutModLengthNorth = (256.54, 309.9, 309.9, 508.19, 508.19, 508.19) #6 rows
MINOSNORTHY = -0.5*SHELLY+0.5*(len(minosCutModLengthNorth)*mModW+(len(minosCutModLengthNorth)-1)*SIDECRTSHELFTHICK+2*PADTagger)-PADTagger+SIDECRTNORTHWALLTOFLOOR-BOTTOMCRTROLLERHEIGHT*0.9

########################################################

gdml = ET.Element('gdml')
if testmode: materials = ET.SubElement(gdml, 'materials')
solids = ET.SubElement(gdml, 'solids')
structure = ET.SubElement(gdml, 'structure')
solids_store = {}
modToFeb = dict()
feb_id = 0
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

def strip(style="m", modnum=0, stripnum=0, length=0):
    '''Build one scintillator strip.'''

    if style=="m":
        x=XM
        y=YM
        if length==0:
            z=ZM
        else:
            z=length
        name = 'MINOS'
    if style=="c":
        name = 'CERN'
        x=XC
        if stripnum < NXC:
            y=YCTOP
        else:
            y=YCBOT
        z=ZC
    if style=="d":
        x=XD
        y=YD
        z=ZD
        name = 'DC'

    xx = str(x)
    yy = str(y)
    zz = str(z) 

    sname = 'AuxDetSensitive_' + name
    vname = 'volAuxDetSensitive_' + name + '_module_'

    if modnum < 10:
        vname += '00'
    elif modnum < 100: 
        vname += '0'

    vname += str(modnum)+'_'
    if style=='m' and length!=0:
        sname+='_cut'+str(int(length))+'_'
        vname+='cut'+str(int(length))+'_'

    if style=='c':
        if stripnum < NXC:
            sname += '_top_'
            vname += 'top_'
        else:
            sname += '_bot_'
            vname += 'bot_'

    vname+='strip_'
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

def module(style="c", reg='tt', length=0):
    '''Build an edge-to-edge array of scintillator strips.'''

    if style=="m":
        ny=NXM
        x=XM
        y=YM
        if length==0:
            z=ZM
            zz = str(mModL)
            zzsub = str(mModL-2*PADM)
        else:
            if length>ZM:
                print('WARNING: LENGTH SPECIFIED FOR CUT MINOS MODULE EXCEEDES FULL LENGTH!')
            z  = length
            zz = str(mModL-ZM+length)
            zzsub = str(mModL-2*PADM-ZM+length)
        xx = str(mModH) 
        yy = str(mModW) 
        xxsub = str(mModH-2*PADM) 
        yysub = str(mModW-2*PADM) 
        name = "MINOS"

    if style=="c":
        x=XC
        z=ZC
        ny=NXC
        xx = str(cModW) 
        yy = str(cModH) 
        zz = xx
        xxsub = str(cModW - 2*PADC) 
        yysub = str(cModH-2*PADC)
        zzsub = xxsub
        name = "CERN"

    if style=="d":
        x=XD
        y=YD
        z=ZD
        ny=NXD
        xx = str(x*(ny+0.5)+2*PADD+(ny+2)*PADStrip)
        yy = str(dModH) 
        zz = str(z+2*PADD+2*PADStrip)
        xxsub = str(x*(ny+0.5)+(ny+2)*PADStrip)
        yysub = str(dModH-2*PADM) 
        zzsub = str(z+2*PADStrip)
        name = "DC"

    modnum = get_mod_id(style)
    stripnum = 0

    sname = 'AuxDet_' + name + '_module_'
    vname = 'vol' + sname

    if modnum < 10:
        vname += '00'
    elif modnum < 100:
        vname += '0'
    vname += str(modnum)+'_'

    if style=='m' and length!=0:
        sname += 'cut'+str(int(length))
        vname += 'cut'+str(int(length))+'_'

    if reg=='tt': vname += 'Top'
    if reg=='rn': vname += 'RimNorth' 
    if reg=='rs': vname += 'RimSouth' 
    if reg=='rw': vname += 'RimWest' 
    if reg=='re': vname += 'RimEast' 
    if reg=='ss': vname += 'South' 
    if reg=='nn': vname += 'North' 
    if reg=='ws': vname += 'WestSouth' 
    if reg=='wc': vname += 'WestCenter'
    if reg=='wn': vname += 'WestNorth'        
    if reg=='es': vname += 'EastSouth' 
    if reg=='ec': vname += 'EastCenter'
    if reg=='en': vname += 'EastNorth'
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

    #generate strips, top layer for c or d type
    for i in range(ny):
        if style=='m' and length>0:
            strips.append(strip(style, modnum, stripnum,length))
        else:
            strips.append(strip(style, modnum, stripnum,0))
        stripnum += 1

    #generate bottom layer strips
    if style=='d' or style=='c':    
        for i in range(ny):
            strips2.append(strip(style, modnum, stripnum))
            stripnum += 1

    vnamein = vname + '_inner'
    vin = ET.SubElement(structure, 'volume', name=vnamein)
    ET.SubElement(vin, 'materialref', ref='Air')
    ET.SubElement(vin, 'solidref', ref=snamein)

    #place first layer of strips (only layer for m modules)
    #top layer for c or d modules
    for i, (es, ev) in enumerate(strips):
        pv = ET.SubElement(vin, 'physvol')
        ET.SubElement(pv, 'volumeref', ref=ev.attrib['name'])

        if style=='m':
            dy = (2*i - ny + 1)* 0.5 * (y+PADStrip)
            dx=0
        if style=='c':
            dx = (2*i - ny + 1)* 0.5 * (x+PADStrip)
            dy=0.5*(YCBOT+PADStrip)

        if style=='d':
            dy= 0.5*(y+PADStrip)
            dx=(i - 0.5*ny + 0.25) * (x+PADStrip)

        posname = 'pos' + ev.attrib['name']
        ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(dx), y=str(dy), z='0')

    #place bottom layers
    if style=='c':
        for i, (es, ev) in enumerate(strips2):
            pv = ET.SubElement(vin, 'physvol')
            ET.SubElement(pv, 'volumeref', ref=ev.attrib['name'])

            dy= -0.5*(YCTOP+PADStrip)
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

#build one stack of MINOS modules for east or west side CRT walls
#pos specifies one of 3 stacks on either side: 's, 'c', 'n'
def minosSideTagger(side='e', pos='n'): 
    ''' Build a side tagger (1 stacks)
    '''
    if side!='e' and side !='w':
        print('bad value passed to function, minosSideTagger: side='+side)
    if pos!='n' and pos!='c' and pos != 's':
        print('bad value passed to function, minosSideTagger: pos='+pos)

    coords = []
    modules = []
    nstack=NMODSTACK

    if(pos=='c'):
        nstack-=1

    z = mModL+2*PADTagger
    if(pos=='s'):
        z+=MINOSLATSOUTHACTIVEOVERHANG

    xx = str(SIDECRTSTACKDEPTH)
    yy = str(nstack*mModW+(nstack-1)*SIDECRTSHELFTHICK+2*PADTagger)
    zz = str(z)

    if pos=='s':
        xxsub = str(SIDECRTSTACKDEPTH)
        yysub = str(mModW+PADTagger)
        zzsub = str(MINOSLATSOUTHACTIVEOVERHANG+PADTagger)
        xpsub = '0'
        ypsub = str(0.5*((nstack-1)*(mModW+SIDECRTSHELFTHICK)+PADTagger))
        zpsub = str(-0.5*(z-MINOSLATSOUTHACTIVEOVERHANG))

    #loop over stacks
    for layer in range (2): #6):
        
        dx = ((-1)**layer)*0.5*LAYERSPACE
        if side=='w': dx*= -1
        dz=0
        #loop over modules in stack
        for i in range(nstack):
            dy = 0.5*(2*i+1-nstack)*(mModW+SIDECRTSHELFTHICK)
            #if pos=='c':
            #    dy-=0.5*mModW
            if pos=='s' and i==nstack-1 : dz = 0.5*MINOSLATSOUTHACTIVEOVERHANG
            elif pos=='s': dz = -0.5*MINOSLATSOUTHACTIVEOVERHANG
            coords.append((dx,dy,dz))

    sname = 'tagger_SideLat_'
    if pos=='c': sname+='Center'
    if pos=='s': 
        sname+='South'
        snameint = sname + '_internal'
        snameext = sname + '_external'
    if pos=='n': sname+='North'

    if not sname in solids_store:
        if pos != 's': stagger = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
        if pos=='s':
            sext = ET.SubElement(solids, 'box', name=snameext, lunit="cm", x=xx, y=yy, z=zz)
            sint = ET.SubElement(solids, 'box', name=snameint, lunit="cm", x=xxsub, y=yysub, z=zzsub)
            stagger = ET.SubElement(solids, 'subtraction', name=sname)
            ET.SubElement(stagger, 'first', ref=snameext)
            ET.SubElement(stagger, 'second', ref=snameint)
            ET.SubElement(stagger, 'position', name='crtsouthtaggersubpos', unit='cm', x=xpsub,y=ypsub,z=zpsub)
            solids_store[snameext] = sext
            solids_store[snameint] = sint
        solids_store[sname] = stagger
    else:
        if pos=='s':
            sext = solids_store[snameext]
            sint = solids_store[snameint]
        stagger = solids_store[sname]

    vname = 'vol_'+ sname+'_'

    global feb_id
    global printModIds
    fmod = 0

    if printModIds: print('MINOS tagger, '+side+', pos '+pos+' first module: '+str(mod_id+1)+', first FEB: '+str(feb_id+1))
    feb_id+=2

    for i in range(len(coords)):
        if side=='w':
            if pos=='s':
                modules.append(module('m','ws'))
                if i==0: vname+='WestSouth'
            if pos=='c':
                modules.append(module('m','wc'))
                if i==0: vname+='WestCenter'
            if pos=='n':
                modules.append(module('m','wn'))
                if i==0: vname+='WestNorth'
        if side=='e':
            if pos=='s':
                modules.append(module('m','es'))
                if i==0: vname+='EastSouth'
            if pos=='c':
                modules.append(module('m','ec'))
                if i==0: vname+='EastCenter'
            if pos=='n':
                modules.append(module('m','en'))
                if i==0: vname+='EastNorth'
        fmod+=1
        modToFeb[mod_id] = ((feb_id-1,fmod),(feb_id,fmod))
        if fmod==3:
            fmod=0
            if i != len(coords)-1: feb_id+=2

    if printModIds: print('   last module: '+str(mod_id)+', last FEB: '+str(feb_id))
    #print('dictionary generated:')
    #for k in modToFeb.keys():
    #    print('module '+str(k)+': '+str(modToFeb[k]))

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

def minosNorthTagger():
    ''' Build north MINOS tagger (2 layers cut modules in X-X) on upstream face
    '''

    coords = []
    modules = []
    ny = len(minosCutModLengthNorth)

    x  = 2*(mModL-ZM+max(minosCutModLengthNorth))+PADM+2*PADTagger
    y  = ny*mModW+(ny-1)*SIDECRTSHELFTHICK+2*PADTagger
    xx = str(x)
    yy = str(y)
    zz = str(SIDECRTSTACKDEPTH)

    global feb_id
    global printModIds
    if printModIds: print('MINOS tagger North, first module: '+str(mod_id+1)+', FEB: '+str(feb_id+1))
    fmod = 0
    feb_id+=4

    #loop over rows starting from bottom going up
    for row in range(ny):
        zin   = 0.5*LAYERSPACE
        xleft = 0.5*x - PADTagger - 0.5*(mModL-ZM+minosCutModLengthNorth[row])
        yrow  = -0.5*y + PADTagger + (row+0.5)*mModW + row*SIDECRTSHELFTHICK
        rowcoords = ( (xleft,yrow,zin),(-xleft,yrow,zin),(xleft,yrow,-zin),(-xleft,yrow,-zin) )
        coords.append(rowcoords)
        rowmodules = []
        fmod+=4
        for i in range(4):
            rowmodules.append(module('m','nn',minosCutModLengthNorth[row]))
            modToFeb[mod_id] = (feb_id-3+i,fmod/4)
        modules.append(rowmodules)
        if fmod==12:
            fmod=0
            if row != ny-1: feb_id+=4

    if printModIds: print('   last module: '+str(mod_id)+', FEB: '+str(feb_id))

    sname = 'tagger_SideNorth'
    vname = 'vol_'+ sname

    stagger = ET.SubElement(solids, 'box', name=sname, lunit="cm", x=xx, y=yy, z=zz)
    vtagger = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(vtagger, 'materialref', ref='Air')
    ET.SubElement(vtagger, 'solidref', ref=sname)

    #place left side module phy. vol.s
    for row in range(len(coords)):
        for mod, (xc,yc,zc) in enumerate(coords[row]):

            (s,v)=modules[row][mod]
            pv = ET.SubElement(vtagger, 'physvol')
            ET.SubElement(pv, 'volumeref', ref=v.attrib['name'])

            posname = 'pos' + v.attrib['name']
            ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

            if xc>0: 
                posname = 'rotplus' + v.attrib['name']
                ET.SubElement(pv, 'rotation', name=posname, unit="deg", x='0', y='90', z='0')
            else:
                posname = 'rotneg' + v.attrib['name'] 
                ET.SubElement(pv, 'rotation', name=posname, unit="deg", x='0', y='-90', z='0')

    return stagger, vtagger

def minosSouthTagger(): 
    ''' Build front MINOS tagger (2 layers in X-Y) on upstream face
    '''
    nmody = 9
    coords = []
    modules = []

    x = mModL+SIDECRTSOUTHWALLLATOFFSET+2*PADTagger
    y = 9*mModW+8*SIDECRTSHELFTHICK+2*PADTagger
    z = SIDECRTSTACKDEPTH+mModH+PADTagger

    xx = str(x)
    yy = str(y)
    zz = str(z)

    for i in range(2*nmody):

        if i < nmody: #bottom row
            dx = 0.5*x-PADTagger-SIDECRTSOUTHWALLLATOFFSET-(i+0.5)*mModW-i*PADModule
            dy = -0.5*y + PADTagger + 0.5*(mModL-0.5*ZM)
            dz = -0.5*z + PADTagger + 0.5*mModH
        else: #top row
            dx = 0.5*x-PADTagger-SIDECRTSOUTHWALLLATOFFSET-(i+0.5-nmody)*mModW-(i-nmody)*PADModule
            dy = 0.5*y - PADTagger - 0.5*(mModL - 0.5*ZM)
            dz = -0.5*z +PADTagger + mModH + 1.5*SIDECRTPOSTWIDTH

        coords.append((dx,dy,dz,1)) #x,y,z,vert=true


    for i in range(NMODSTACK):

        dx = 0.5*x-PADTagger-0.5*mModL
        dy = -0.5*y+PADTagger+(i+0.5)*mModW + i*SIDECRTSHELFTHICK
        if i > NMODSTACK-3:
            dx -= SIDECRTSOUTHWALLLATOFFSET 
        dz = 0.5*z - 1.5*SIDECRTPOSTWIDTH

        coords.append((dx,dy,dz,0)) #x,y,z,vert=false

    global feb_id
    global printModIds
    if printModIds: print('MINOS tagger South, first module: '+str(mod_id+1)+', FEB: '+str(feb_id+1))
    fmod = 0
    feb_id+=1

    for i in range(len(coords)):
        if i<2*nmody:
            modules.append(module('m','ss',0.5*ZM))
            fmod+=1
            modToFeb[mod_id] = (feb_id,fmod)
            if fmod==3:
                fmod=0
                if i!=2*nmody-1: feb_id+=1
                else: feb_id+=2
        else:
            modules.append(module('m','ss',0))
            fmod+=1
            modToFeb[mod_id] = ((feb_id-1,fmod),(feb_id,fmod))
            if fmod==3:
                fmod=0
                if i!= len(coords)-1: feb_id+=2

    if printModIds: print('   last module: '+str(mod_id)+', FEB: '+str(feb_id))

    sname = 'tagger_SideSouth'
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


def DCTagger():
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

        coords.append((dx,0,dz,rot))

    global feb_id
    global printModIds
    if printModIds: print('DC tagger, first module: '+str(mod_id+1)+', FEB: '+str(feb_id+1))

    for i in range(len(coords)):
        modules.append(module('d','bt'))
        feb_id+=1
        modToFeb[mod_id] = (feb_id,1)

    if printModIds: print('   last module: '+str(mod_id)+', FEB: '+str(feb_id))

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


def cernTopTagger():
    ''' Build top CERN tagger (1 layer of modules) 
    '''
    modwidth = ZC + 2*PADC + (NXC+1)*PADStrip
    xx = str(NTOPX*modwidth+2*PADTagger+(NTOPX-1)*PADModule)
    yy = str(YCTOP+YCBOT+3*PADStrip+2*PADC+2*PADTagger)
    zz = str(NTOPZ*modwidth + 2*PADTagger + (NTOPZ-1)*PADModule)

    coords = []
    modules = []

    dz = 0.5*(modwidth+PADModule)*(1 - NTOPZ)
    dx = 0.5*(modwidth+PADModule)*(1 - NTOPX)

    for i in range(NTOPX*NTOPZ):

        coords.append((dx,0,dz))

        if (i+1)%NTOPZ == 0:
            dx+= modwidth + PADModule
            dz = 0.5*(modwidth+PADModule)*(1 - NTOPZ)
        else: dz+= modwidth + PADModule

    global feb_id
    global printModIds
    if printModIds: print('CERN tagger Top, first module: '+str(mod_id+1)+', FEB: '+str(feb_id+1))

    for i in range(len(coords)):
        modules.append(module('c','tt'))
        feb_id+=1
        modToFeb[mod_id] = (feb_id,1)

    if printModIds: print('   last module: '+str(mod_id)+', FEB: '+str(feb_id))

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

def cernLatRimTagger(side='L'):
    ''' Build east(side='R') or west(side='L') CERN rim tagger
    '''
    modwidth = ZC + 2*PADC + (NXC+1)*PADStrip
    xx = str(modwidth+2*PADTagger)
    yy = str(YCTOP+YCBOT+3*PADStrip+2*PADC+2*PADTagger)
    zz = str(NSLOPELAT*modwidth + 2*PADTagger + (NSLOPELAT-1)*PADModule)

    coords = []
    modules = []

    dz = 0.5*(modwidth+PADModule)*(1 - NSLOPELAT)

    for i in range(NSLOPELAT):

        coords.append((0,0,dz))

        dz+= modwidth+PADModule

    global feb_id
    global printModIds
    if printModIds: print('CERN tagger Lat, side '+side+' first module: '+str(mod_id+1)+', FEB: '+str(feb_id+1))
    
    for i in range(len(coords)):
        if side == 'L':
            modules.append(module('c','rw'))
        if side == 'R':
            modules.append(module('c','re'))
        feb_id+=1
        modToFeb[mod_id] = (feb_id,1)

    if printModIds: print('   last module: '+str(mod_id)+', FEB: '+str(feb_id))

    sname = 'tagger_'
    if side == 'L':
        sname += 'RimWest'
    if side == 'R':
        sname += 'RimEast'
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


def cernLongRimTagger(side='U'):
    ''' Build  north(side='D') or south(side='U') CERN rim tagger
    '''
    modwidth = ZC + 2*PADC + (NXC+1)*PADStrip
    xx = str(NSLOPEFRONT*modwidth+2*PADTagger+(NSLOPEFRONT-1)*PADModule)
    yy = str(YCTOP+YCBOT+3*PADStrip+2*PADC+2*PADTagger)
    zz = str(modwidth + 2*PADTagger)

    coords = []
    modules = []

    dx = 0.5*(modwidth+PADModule)*(1 - NSLOPEFRONT)

    for i in range(NSLOPEFRONT):

        coords.append((dx,0,0))
        dx+= modwidth+PADModule

    global feb_id
    global printModIds
    if printModIds: print('CERN tagger Long, side '+side+' first module: '+str(mod_id+1)+', FEB: '+str(feb_id+1))

    for i in range(len(coords)):
        if side == 'U':
            modules.append(module('c','rs'))
        if side == 'D':
            modules.append(module('c','rn'))
        feb_id+=1
        modToFeb[mod_id] = (feb_id,1)

    if printModIds: print('   last module: '+str(mod_id)+', FEB: '+str(feb_id))

    sname = 'tagger_'
    if side == 'U':
        sname += 'RimSouth'
    if side == 'D':
        sname += 'RimNorth'
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

    #shell outer and void dimensions
    WVPADY = 25
    xxint = str(WVWIDTH + 2*SIDECRTWVOFFSET)
    yyint = str(WVHEIGHT+1.0+WVPADY) 
    zzint = str(WVLENGTH+2*SIDECRTWVOFFSET) 

    xxext = str(WVWIDTH + 2*SIDECRTROLLOFFSET + 1.1*SIDECRTSTACKDEPTH) 
    yyext = str(SHELLY) 
    zzext = str(SHELLZ) 

    #generate all of the tagger volumes, CRT modules, and strips
    (s,vws) = minosSideTagger('w','s') #MINOS west wall, south stack
    (s,vwc) = minosSideTagger('w','c') #MINOS west wall, center stack
    (s,vwn) = minosSideTagger('w','n') #MINOS west wall, north stack
    (s,ves) = minosSideTagger('e','s') #MINOS east wall, south stack
    (s,vec) = minosSideTagger('e','c') #MINOS east wall, center stack
    (s,ven) = minosSideTagger('e','n') #MINOS east wall, north stack
    (s,vss) = minosSouthTagger()#'U',0,0,0) #MINOS south
    (s,vnn) = minosNorthTagger() #MINOS north
    (s,vbt) = DCTagger() #DC Bottom
    (s,vtt) = cernTopTagger() #CERN top
    (s,vrw) = cernLatRimTagger('L') #CERN RimWest
    (s,vre) = cernLatRimTagger('R') #CERN RimEast
    (s,vrs) = cernLongRimTagger('U') #CERN RimSouth
    (s,vrn) = cernLongRimTagger('D') #CERN RimNorth

    #CRT Shell containing all of the tagger volumes and a void to cointain the warm vessel
    sname = 'CRT_Shell'
    snameext = sname+'_external'
    snameint = sname+'_internal'
    sexternal = ET.SubElement(solids, 'box', name=snameext, lunit="cm", x=xxext, y=yyext, z=zzext)
    sinternal = ET.SubElement(solids, 'box', name=snameint, lunit="cm", x=xxint, y=yyint, z=zzint)
    sshell = ET.SubElement(solids, 'subtraction', name=sname)
    ET.SubElement(sshell, 'first', ref=snameext)
    ET.SubElement(sshell, 'second', ref=snameint)
    ET.SubElement(sshell, 'position', name='crtshellsubpos', unit='cm', x='0',y=str(-0.5*SHELLY+WVFOOTELEVATION+0.5*WVHEIGHT+0.5*WVPADY),z=str(-SHELLWVOFFSET))

    vname = 'vol'+sname
    vshell = ET.SubElement(structure, 'volume', name=vname)
    ET.SubElement(vshell, 'materialref', ref='Air')
    ET.SubElement(vshell, 'solidref', ref=sname)

    #position MINOS west, south
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vws.attrib['name'])

    xc = 0.5* WVWIDTH + SIDECRTWVOFFSET + 0.5*SIDECRTSTACKDEPTH 
    yc = MINOSLATFIXY 
    zc = MINOSLATSOUTHZ - SHELLWVOFFSET 

    posname = 'pos' + vws.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    #position MINOS west, center
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vwc.attrib['name'])

    xc = 0.5* WVWIDTH + SIDECRTROLLOFFSET 
    yc = MINOSLATROLLY 
    zc = MINOSLATCENTZ - SHELLWVOFFSET 

    posname = 'pos' + vwc.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    #position MINOS west, north
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vwn.attrib['name'])

    xc = 0.5* WVWIDTH + SIDECRTWVOFFSET + 0.5*SIDECRTSTACKDEPTH 
    yc = MINOSLATFIXY 
    zc = MINOSLATNORTHZ - SHELLWVOFFSET 

    posname = 'pos' + vwn.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    #Position MINOS east, south
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=ves.attrib['name'])

    xc = -0.5* WVWIDTH - SIDECRTWVOFFSET - 0.5*SIDECRTSTACKDEPTH 
    yc = MINOSLATFIXY 
    zc = MINOSLATSOUTHZ - SHELLWVOFFSET 

    posname = 'pos' + ves.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    #Position MINOS east, center
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vec.attrib['name'])

    xc = -0.5* WVWIDTH - SIDECRTROLLOFFSET 
    yc = MINOSLATROLLY 
    zc = MINOSLATCENTZ - SHELLWVOFFSET 

    posname = 'pos' + vec.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    #Position MINOS east, north
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=ven.attrib['name'])

    xc = -0.5* WVWIDTH - SIDECRTWVOFFSET - 0.5*SIDECRTSTACKDEPTH 
    yc = MINOSLATFIXY 
    zc = MINOSLATNORTHZ - SHELLWVOFFSET 

    posname = 'pos' + ven.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    #position MINOS north (downstream)
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vnn.attrib['name'])

    xc = 0.0 
    yc = MINOSNORTHY 
    zc = 0.5* WVLENGTH + SIDECRTWVOFFSET + 0.5*SIDECRTSTACKDEPTH - SHELLWVOFFSET 

    posname = 'pos' + vnn.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    #position MINOS south (upstream)
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vss.attrib['name'])

    xc = 0.0 
    yc = MINOSSOUTHY 
    zc = -1*(0.5* WVLENGTH + SIDECRTWVOFFSET + 0.5*(SIDECRTSTACKDEPTH+PADTagger+mModH)) - SHELLWVOFFSET 

    posname = 'pos' + vss.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    #position DC Bottom
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vbt.attrib['name'])

    xc = posDCInDetEncl[0]
    yc = posDCInDetEncl[1]
    zc = posDCInDetEncl[2] - SHELLWVOFFSET

    posname = 'pos' + vbt.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))  

    #position CERN Top
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vtt.attrib['name'])

    xc = 0.0
    yc = CERNTOPY
    zc = CERNTOPZ - SHELLWVOFFSET 

    posname = 'pos' + vtt.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc)) 
   
    #position CERN west rim
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vrw.attrib['name'])

    xc = CERNRIMLATX 
    yc = CERNRIMLATY 
    zc = CERNRIMLATZ - SHELLWVOFFSET 

    posname = 'pos' + vrw.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc)) 

    posname = 'rot' + vrw.attrib['name']
    ET.SubElement(pv, 'rotation', name=posname, unit="deg", x='0', y='0', z=str(SLOPEINCLINATION))
 
    #position CERN east rim
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vre.attrib['name'])

    xc = -1*CERNRIMLATX 
    yc = CERNRIMLATY 
    zc = CERNRIMLATZ - SHELLWVOFFSET 

    posname = 'pos' + vre.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))

    posname = 'rot' + vre.attrib['name']
    ET.SubElement(pv, 'rotation', name=posname, unit="deg", x='0', y='0', z=str(-1*SLOPEINCLINATION))

    #position CERN south rim
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vrs.attrib['name'])

    xc = 0.0 
    yc = CERNRIMSY 
    zc = CERNRIMSZ - SHELLWVOFFSET 

    posname = 'pos' + vrs.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))
    
    posname = 'rot' + vrs.attrib['name']
    ET.SubElement(pv, 'rotation', name=posname, unit="deg", x=str(SLOPEINCLINATION), y='0', z='0')

    #position CERN north rim
    pv = ET.SubElement(vshell, 'physvol')
    ET.SubElement(pv, 'volumeref', ref=vrn.attrib['name'])

    xc = 0.0 
    yc = CERNRIMNY 
    zc = CERNRIMNZ - SHELLWVOFFSET 

    posname = 'pos' + vrn.attrib['name']
    ET.SubElement(pv, 'position', name=posname, unit="cm", x=str(xc), y=str(yc), z=str(zc))
   
    posname = 'rot' + vrn.attrib['name']
    ET.SubElement(pv, 'rotation', name=posname, unit="deg", x=str(-1*SLOPEINCLINATION), y='0', z='0')

    return sshell,vshell

#######################################################################
# now if test mode generate materials, CRT shell, world, gdml header
# else just generate CRT shell for use with master geometry generator script

if testmode:
   m = ET.SubElement(materials, 'element', name='aluminum', formula='Al', Z='13')
   ET.SubElement(m, 'atom', value='26.9815')

   m = ET.SubElement(materials, 'element', name='nitrogen', formula='N', Z='7')
   ET.SubElement(m, 'atom', value='14.0067')

   m = ET.SubElement(materials, 'element', name='oxygen', formula='O', Z='8')
   ET.SubElement(m, 'atom', value='15.999')

   m = ET.SubElement(materials, 'element', name='argon', formula='Ar', Z='18')
   ET.SubElement(m, 'atom', value='39.9480')

   m = ET.SubElement(materials, 'element', name='hydrogen', formula='H', Z='1')
   ET.SubElement(m, 'atom', value='1.0079')

   m = ET.SubElement(materials, 'element', name='carbon', formula='C', Z='6')
   ET.SubElement(m, 'atom', value='12.0107')

   m = ET.SubElement(materials, 'material', name='ALUMINUM_Al', formula='ALUMINUM_Al')
   ET.SubElement(m, 'D', value='2.6990', unit='g/cm3')
   ET.SubElement(m, 'fraction', n='1.000', ref='aluminum')

   m = ET.SubElement(materials, 'material', name='Air')
   ET.SubElement(m, 'D', value='0.001205', unit='g/cm3')
   ET.SubElement(m, 'fraction', n='0.781154', ref='nitrogen')
   ET.SubElement(m, 'fraction', n='0.209476', ref='oxygen')
   ET.SubElement(m, 'fraction', n='0.00934', ref='argon') 

   m  = ET.SubElement(materials, 'material', name='Polystyrene')
   ET.SubElement(m, 'D', value='1.19', unit='g/cm3')
   ET.SubElement(m, 'fraction', n='0.077418', ref='hydrogen')
   ET.SubElement(m, 'fraction', n='0.922582', ref='carbon')

#crt shell volume
(s,v) = detectorEnclosure()
#(s,v) = module()

if testmode:
   ws = ET.SubElement(solids, 'box', name='World', lunit="cm", x='1500', y='1500', z='3000')
   w = ET.SubElement(structure, 'volume', name='volWorld')
   ET.SubElement(w, 'materialref', ref='Air')
   ET.SubElement(w, 'solidref', ref='World')
   pv = ET.SubElement(w, 'physvol')
   ET.SubElement(pv, 'volumeref', ref=v.attrib['name'])
   posname = 'pos' + v.attrib['name']
   ET.SubElement(pv, 'position', name=posname, unit="cm", x='0', y='0', z='0')
   setup = ET.SubElement(gdml, 'setup', name='Default', version='1.0')
   ET.SubElement(setup, 'world', ref='volWorld')

#no. modules generated for each subsystem
print('MINOS modules generated: '+str(nModM))
print('CERN  modules generated: '+str(nModC))
print('DblCh modules generated: '+str(nModD))

#write to file
with open('icarus_crt.gdml', 'w') as f:
    f.write(minidom.parseString(ET.tostring(gdml)).toprettyxml(indent='\t'))

print('Writing dictionary to file...')
with open('feb_map.txt','w') as f:
    writer = csv.writer(f,delimiter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
    for mod in modToFeb.keys():
        tmp = modToFeb[mod]
        if type(tmp[0])==int:
            writer.writerow([str(mod),str(tmp[0]),str(tmp[1])])
        else:
             writer.writerow([str(mod),str(tmp[0][0]),str(tmp[0][1]),str(tmp[1][0]),str(tmp[1][1])])
