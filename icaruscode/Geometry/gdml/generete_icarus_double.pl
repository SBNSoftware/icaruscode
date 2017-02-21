#!/usr/bin/perl

# contact tylerdalion@gmail.com for any GDML/generate questions
# I would love to help!


#
# Basic ICARUS detector, much tweaking to do.
#


# Each subroutine generates a fragment GDML file, and the last subroutine
# creates an XML file that make_gdml.pl will use to appropriately arrange
# the fragment GDML files to create the final desired DUNE GDML file, 
# to be named by make_gdml output command

# If you are playing with different geometries, you can use the
# suffix command to help organize your work.

use Math::Trig;
use Getopt::Long;
use Math::BigFloat;
Math::BigFloat->precision(-16);

GetOptions( "help|h" => \$help,
	    "suffix|s:s" => \$suffix,
	    "output|o:s" => \$output,
	    "wires|w:s" => \$wires,
            "helpcube|c" => \$helpcube);

if ( defined $help )
{
    # If the user requested help, print the usage notes and exit.
    usage();
    exit;
}

if ( ! defined $suffix )
{
    # The user didn't supply a suffix, so append nothing to the file
    # names.
    $suffix = "";
}
else
{
    # Otherwise, stick a "-" before the suffix, so that a suffix of
    # "test" applied to filename.gdml becomes "filename-test.gdml".
    $suffix = "-" . $suffix;
}



#++++++++++++++++++++++++ Begin defining variables +++++++++++++++++++++++++

# Define detector geometry variables - later to be put in a parameters
# XML file to be parsed as an input?

# set wires on to be the default, unless given an input by the user
$wires_on = 1; # 1=on, 0=off
if (defined $wires)
{
    $wires_on = $wires;
}

$tpc_on=1;
$inch = 2.54;


##################################################################
##################### wire plane parameters ######################

$UWirePitch             =   0.3;
$VWirePitch             =   0.3;
$XWirePitch             =   0.3;

$UAngle                 =   60;
$VAngle                 =   60;

#prova veloce
#$UAngle                 =   30;
#$VAngle                 =   30;


$SinUAngle              =   sin( deg2rad($UAngle) );
$CosUAngle              =   cos( deg2rad($UAngle) );
$TanUAngle              =   tan( deg2rad($UAngle) );

$SinVAngle              =   sin( deg2rad($VAngle) );
$CosVAngle              =   cos( deg2rad($VAngle) );
$TanVAngle              =   tan( deg2rad($VAngle) );

$UWireCornerInt_y       =   $UWirePitch * $CosUAngle;
$UWireCornerInt_z       =   $UWirePitch * $SinUAngle;
$UWire_ypitch           =   $UWirePitch / $SinUAngle;
$UWire_zpitch           =   $UWirePitch / $CosUAngle;


$VWireCornerInt_y       =   $VWirePitch * $CosVAngle;
$VWireCornerInt_z       =   $VWirePitch * $SinVAngle;
$VWire_ypitch           =   $VWirePitch / $SinVAngle;
$VWire_zpitch           =   $VWirePitch / $CosVAngle;


$TPCWireThickness       =   0.015;

$TPCWirePlane_x  =   $TPCWireThickness;
#height and length defined lower


$nSkipCornerUWires = 0;
$nSkipCornerVWires = 0;


###########################################################################
########################### spacing parameters ############################


$CPA_x                 =     0.51; 
$WirePlaneSpacing      =     0.3;   # center to center
 # NOTE: if wire plane spacing is change, GeometeryTest must be modified
$MaxDrift              =     150; 
 #MaxDrift is the distance form the edge of the CPA to the edge of the first wire plane
 
#Cryostat space with LAr outside of entire fiducial volume
$SpaceWirePlToWall     =     50; 
$SpaceWirePlToWirePl   =     85; # edge to edge, like all of these
$SpaceTPCToFloor       =     50; 
$SpaceTPCToTopLAr      =     50;  
$UpstreamLArPadding    =     50;
$DownstreamLArPadding  =     50;



##############################################################
############## Cryo and TPC relevant dimensions  #############


$TPC_x    =     $MaxDrift + $TPCWirePlane_x + 3*$WirePlaneSpacing;
$TPC_y    =     390;
$TPC_z    =     1990;


$LAr_x    =     $CPA_x 
              + 2*($TPC_x + $SpaceWirePlToWall);
$LAr_y    =     $TPC_y 
              + $SpaceTPCToFloor 
              + $SpaceTPCToTopLAr;
$LAr_z    =     $TPC_z
              + $UpstreamLArPadding 
              + $DownstreamLArPadding;

$SteelThickness		=	0.5*$inch; #half inch
$GaseousAr_y            =       20;

$Cryostat_x = $LAr_x + 2*$SteelThickness ;
$Cryostat_y = $LAr_y + 2*$SteelThickness + $GaseousAr_y ;
$Cryostat_z = $LAr_z + 2*$SteelThickness ;

$LAr_x_orig =   4*($TPC_x)                                         #for total positioning
              + 2*($CPA_x + $SpaceWirePlToWall) 
              + $SpaceWirePlToWirePl;

$Cryostat_x_orig = $LAr_x_orig + 2*$SteelThickness ;               #for total positioning


# $TPCWirePlane_x  defined higher up
$TPCWirePlane_y     =       $TPC_y ; 
$TPCWirePlane_z     =       $TPC_z ;


$TPCinCryo_x[0]     =      - $TPC_x/2 - $CPA_x/2;
#$TPCinCryo_x[1]     =      - $LAr_x/2 + $SpaceWirePlToWall + 1.5*($TPC_x) + $CPA_x ;
#$TPCinCryo_x[2]     =        $LAr_x/2 - $SpaceWirePlToWall - 1.5*($TPC_x) - $CPA_x ;
$TPCinCryo_x[1]     =        $TPC_x/2 + $CPA_x/2;

$posCat_x      =      0;

#$posRightCat_x      =      - $Cryostat_x/2 + 1.0*($TPC_x) + $SpaceWirePlToWall + $CPA_x/2 ;
#$posLeftCat_x       =        $Cryostat_x/2 - 1.0*($TPC_x) - $SpaceWirePlToWall - $CPA_x/2 ;

$TPCinCryo_y        =      - $Cryostat_y/2 + $TPC_y/2 + $SpaceTPCToFloor;  
$TPCinCryo_z        =      - $Cryostat_z/2 + $TPC_z/2 + $UpstreamLArPadding;  

$TPCActive_x        =       $MaxDrift;
$TPCActive_y        =       $TPCWirePlane_y;
$TPCActive_z        =       $TPCWirePlane_z;


$posTPCActive_x     =       $TPC_x/2-$TPCActive_x/2;
$posTPCActive_y     =       0;
$posTPCActive_z     =       0;




##################################################################
############## DetEnc and World relevant parameters  #############

#The padding is the thermal insulation, which is one volume for both T300 modules
$ConcretePadding        =	50;
$FoamPadding            =       80;
$TotalPadding	        =	$ConcretePadding+$FoamPadding;
$DetEnc_x	        =	$Cryostat_x_orig+2*$TotalPadding;
$DetEnc_y	        =	$Cryostat_y+$ConcretePadding; 
                                    # no foam on bottom or top, no concrete on top
$DetEnc_z               =       $Cryostat_z+2*$TotalPadding;


$Cryo1InDetEnc_x     =       -$Cryostat_x/2 - 5 ;
$Cryo2InDetEnc_x     =        $Cryostat_x/2 + 5 ;
$CryoInDetEnc_y     =       -$DetEnc_y/2 + $ConcretePadding + $Cryostat_y/2;
$CryoInDetEnc_z     =       0;


  # We want the world origin to be at the very front of the fiducial volume.
  # move it to the front of the enclosure, then back it up through the concrete/foam, 
  # then through the Cryostat shell, then through the upstream dead LAr (including the
  # dead LAr on the edge of the TPC, but this is covered in $UpstreamLArPadding).
  # This is to be added to the z position of every volume in volWorld

$OriginZSet             =       $DetEnc_z/2 
                              - $TotalPadding 
                              - $SteelThickness 
                              - $UpstreamLArPadding;

$OriginYSet             =       $DetEnc_y/2
                              - $ConcretePadding
                              - $SteelThickness
                              - $SpaceTPCToFloor
                              - $TPC_y/2;

$OriginXSet             =       $LAr_x_orig/2
                              - $SpaceWirePlToWall
                              - 3*$WirePlaneSpacing
                              - $TPCWirePlane_x;


#ORIGINALE
$World_x            =       2*$DetEnc_x;
$World_y            =       2*$DetEnc_y;
$World_z            =       2*$DetEnc_z;

#New->richiesta Christian
#$World_x            =       4*$DetEnc_x;
#$World_y            =       4*$DetEnc_y;
#$World_z            =       4*$DetEnc_z;



#+++++++++++++++++++++++++ End defining variables ++++++++++++++++++++++++++

# run the sub routines that generate the fragments

gen_Define(); 	 # generates definitions at beginning of GDML
gen_Materials(); # generates materials to be used

gen_TPC();	 # generates wires, wire planes, and puts them in volTPC
	         # This is the bulk of the code, and has zero wires option

gen_pmt();       # places pmts

#gen_pmtBig();    #define pmt big

gen_Cryostat();	 # places (2*nAPAWide x nAPAHigh x nAPALong) volTPC,
		 # half rotated 180 about Y
gen_Enclosure(); # places two cryostats and concrete volumes

gen_World();	 # places the enclosure among DUSEL Rock


write_fragments(); # writes the XML input for make_gdml.pl
			# which zips together the final GDML
exit;


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++ usage +++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub usage()
{
    print "Usage: $0 [-h|--help] [-o|--output <fragments-file>] [-s|--suffix <string>]\n";
    print "       if -o is omitted, output goes to STDOUT; <fragments-file> is input to make_gdml.pl\n";
    print "       -s <string> appends the string to the file names; useful for multiple detector versions\n";
    print "       -h prints this message, then quits\n";
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_Define +++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Define()
{

# Create the <define> fragment file name, 
# add file to list of fragments,
# and open it
    $DEF = "icarus_Def" . $suffix . ".gdml";
    push (@gdmlFiles, $DEF);
    $DEF = ">" . $DEF;
    open(DEF) or die("Could not open file $DEF for writing");


print DEF <<EOF;
<?xml version='1.0'?>
<gdml>
<define>

<!--

-->

   <position name="posActiveInTPC"   unit="cm" x="$posTPCActive_x" y="$posTPCActive_y" z="$posTPCActive_z"/>

   <position name="posTPC0inCryo"    unit="cm" x="$TPCinCryo_x[0]" y="$TPCinCryo_y"    z="$TPCinCryo_z" />
   <position name="posCathode"  unit="cm" x="$posCat_x"  y="$TPCinCryo_y"    z="$TPCinCryo_z" />
   <position name="posTPC1inCryo"    unit="cm" x="$TPCinCryo_x[1]" y="$TPCinCryo_y"    z="$TPCinCryo_z" />

   <position name="posCryo1InDetEnc"  unit="cm" x="$Cryo1InDetEnc_x" y="$CryoInDetEnc_y" z="$CryoInDetEnc_z" />
   <position name="posCryo2InDetEnc"  unit="cm" x="$Cryo2InDetEnc_x" y="$CryoInDetEnc_y" z="$CryoInDetEnc_z" />

   <position name="posDetEncInWorld" unit="cm" x="$OriginXSet"     y="$OriginYSet"     z="$OriginZSet"/>


   <position name="posCenter"           unit="cm" x="0" y="0" z="0"/>
   <rotation name="rPlus90AboutX"       unit="deg" x="90" y="0" z="0"/>
   <rotation name="rPlus90AboutY"	unit="deg" x="0" y="90"   z="0"/>
   <rotation name="rPlus90AboutZ"	unit="deg" x="0" y="0"   z="90"/>
   <rotation name="rMinus90AboutY"      unit="deg" x="0" y="270" z="0"/>
   <rotation name="rMinus90AboutYMinus90AboutX"       unit="deg" x="270" y="270" z="0"/>
   <rotation name="rPlusUAngleAboutX"	unit="deg" x="@{[90-$UAngle]}" y="0"   z="0"/>  
   <rotation name="rPlusVAngleAboutX"	unit="deg" x="@{[90+$VAngle]}" y="0"   z="0"/>
   <rotation name="rPlus180AboutY"	unit="deg" x="0" y="180"   z="0"/>
   <rotation name="rIdentity"		unit="deg" x="0" y="0"   z="0"/>
</define>
</gdml>
EOF
    close (DEF);
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++ gen_Materials +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Materials() 
{

# Create the <materials> fragment file name,
# add file to list of output GDML fragments,
# and open it
    $MAT = "icarus_Materials" . $suffix . ".gdml";
    push (@gdmlFiles, $MAT);
    $MAT = ">" . $MAT;
    open(MAT) or die("Could not open file $MAT for writing");


  print MAT <<EOF;
<materials>
  <element name="videRef" formula="VACUUM" Z="1">  <atom value="1"/> </element>
  <element name="bromine" formula="Br" Z="35"> <atom value="79.904"/> </element>
  <element name="hydrogen" formula="H" Z="1">  <atom value="1.0079"/> </element>
  <element name="nitrogen" formula="N" Z="7">  <atom value="14.0067"/> </element>
  <element name="oxygen" formula="O" Z="8">  <atom value="15.999"/> </element>
  <element name="aluminum" formula="Al" Z="13"> <atom value="26.9815"/>  </element>
  <element name="silicon" formula="Si" Z="14"> <atom value="28.0855"/>  </element>
  <element name="carbon" formula="C" Z="6">  <atom value="12.0107"/>  </element>
  <element name="potassium" formula="K" Z="19"> <atom value="39.0983"/>  </element>
  <element name="chromium" formula="Cr" Z="24"> <atom value="51.9961"/>  </element>
  <element name="iron" formula="Fe" Z="26"> <atom value="55.8450"/>  </element>
  <element name="nickel" formula="Ni" Z="28"> <atom value="58.6934"/>  </element>
  <element name="calcium" formula="Ca" Z="20"> <atom value="40.078"/>   </element>
  <element name="magnesium" formula="Mg" Z="12"> <atom value="24.305"/>   </element>
  <element name="sodium" formula="Na" Z="11"> <atom value="22.99"/>    </element>
  <element name="titanium" formula="Ti" Z="22"> <atom value="47.867"/>   </element>
  <element name="argon" formula="Ar" Z="18"> <atom value="39.9480"/>  </element>
  <element name="sulphur" formula="S" Z="16"> <atom value="32.065"/>  </element>
  <element name="phosphorus" formula="P" Z="16"> <atom value="30.973"/>  </element>

  <material name="Vacuum" formula="Vacuum">
   <D value="1.e-25" unit="g/cm3"/>
   <fraction n="1.0" ref="videRef"/>
  </material>

  <material name="ALUMINUM_Al" formula="ALUMINUM_Al">
   <D value="2.6990" unit="g/cm3"/>
   <fraction n="1.0000" ref="aluminum"/>
  </material>

  <material name="SILICON_Si" formula="SILICON_Si">
   <D value="2.3300" unit="g/cm3"/>
   <fraction n="1.0000" ref="silicon"/>
  </material>

  <material name="epoxy_resin" formula="C38H40O6Br4">
   <D value="1.1250" unit="g/cm3"/>
   <composite n="38" ref="carbon"/>
   <composite n="40" ref="hydrogen"/>
   <composite n="6" ref="oxygen"/>
   <composite n="4" ref="bromine"/>
  </material>

  <material name="SiO2" formula="SiO2">
   <D value="2.2" unit="g/cm3"/>
   <composite n="1" ref="silicon"/>
   <composite n="2" ref="oxygen"/>
  </material>

  <material name="Al2O3" formula="Al2O3">
   <D value="3.97" unit="g/cm3"/>
   <composite n="2" ref="aluminum"/>
   <composite n="3" ref="oxygen"/>
  </material>

  <material name="Fe2O3" formula="Fe2O3">
   <D value="5.24" unit="g/cm3"/>
   <composite n="2" ref="iron"/>
   <composite n="3" ref="oxygen"/>
  </material>

  <material name="CaO" formula="CaO">
   <D value="3.35" unit="g/cm3"/>
   <composite n="1" ref="calcium"/>
   <composite n="1" ref="oxygen"/>
  </material>

  <material name="MgO" formula="MgO">
   <D value="3.58" unit="g/cm3"/>
   <composite n="1" ref="magnesium"/>
   <composite n="1" ref="oxygen"/>
  </material>

  <material name="Na2O" formula="Na2O">
   <D value="2.27" unit="g/cm3"/>
   <composite n="2" ref="sodium"/>
   <composite n="1" ref="oxygen"/>
  </material>

  <material name="TiO2" formula="TiO2">
   <D value="4.23" unit="g/cm3"/>
   <composite n="1" ref="titanium"/>
   <composite n="2" ref="oxygen"/>
  </material>

  <material name="FeO" formula="FeO">
   <D value="5.745" unit="g/cm3"/>
   <composite n="1" ref="iron"/>
   <composite n="1" ref="oxygen"/>
  </material>

  <material name="CO2" formula="CO2">
   <D value="1.562" unit="g/cm3"/>
   <composite n="1" ref="iron"/>
   <composite n="2" ref="oxygen"/>
  </material>

  <material name="P2O5" formula="P2O5">
   <D value="1.562" unit="g/cm3"/>
   <composite n="2" ref="phosphorus"/>
   <composite n="5" ref="oxygen"/>
  </material>

  <material formula=" " name="DUSEL_Rock">
    <D value="2.82" unit="g/cm3"/>
    <fraction n="0.5267" ref="SiO2"/>
    <fraction n="0.1174" ref="FeO"/>
    <fraction n="0.1025" ref="Al2O3"/>
    <fraction n="0.0473" ref="MgO"/>
    <fraction n="0.0422" ref="CO2"/>
    <fraction n="0.0382" ref="CaO"/>
    <fraction n="0.0240" ref="carbon"/>
    <fraction n="0.0186" ref="sulphur"/>
    <fraction n="0.0053" ref="Na2O"/>
    <fraction n="0.00070" ref="P2O5"/>
    <fraction n="0.0771" ref="oxygen"/>
  </material> 

  <material name="fibrous_glass">
   <D value="2.74351" unit="g/cm3"/>
   <fraction n="0.600" ref="SiO2"/>
   <fraction n="0.118" ref="Al2O3"/>
   <fraction n="0.001" ref="Fe2O3"/>
   <fraction n="0.224" ref="CaO"/>
   <fraction n="0.034" ref="MgO"/>
   <fraction n="0.010" ref="Na2O"/>
   <fraction n="0.013" ref="TiO2"/>
  </material>

  <material name="FR4">
   <D value="1.98281" unit="g/cm3"/>
   <fraction n="0.47" ref="epoxy_resin"/>
   <fraction n="0.53" ref="fibrous_glass"/>
  </material>

  <material name="STEEL_STAINLESS_Fe7Cr2Ni" formula="STEEL_STAINLESS_Fe7Cr2Ni">
   <D value="7.9300" unit="g/cm3"/>
   <fraction n="0.0010" ref="carbon"/>
   <fraction n="0.1792" ref="chromium"/>
   <fraction n="0.7298" ref="iron"/>
   <fraction n="0.0900" ref="nickel"/>
  </material>

  <material name="LAr" formula="LAr">
   <D value="1.40" unit="g/cm3"/>
   <fraction n="1.0000" ref="argon"/>
  </material>

  <material name="ArGas" formula="ArGas">
   <D value="0.00166" unit="g/cm3"/>
   <fraction n="1.0" ref="argon"/>
  </material>

  <material formula=" " name="Air">
   <D value="0.001205" unit="g/cm3"/>
   <fraction n="0.781154" ref="nitrogen"/>
   <fraction n="0.209476" ref="oxygen"/>
   <fraction n="0.00934" ref="argon"/>
  </material>

  <material formula=" " name="G10">
   <D value="1.7" unit="g/cm3"/>
   <fraction n="0.2805" ref="silicon"/>
   <fraction n="0.3954" ref="oxygen"/>
   <fraction n="0.2990" ref="carbon"/>
   <fraction n="0.0251" ref="hydrogen"/>
  </material>

  <material formula=" " name="Granite">
   <D value="2.7" unit="g/cm3"/>
   <fraction n="0.438" ref="oxygen"/>
   <fraction n="0.257" ref="silicon"/>
   <fraction n="0.222" ref="sodium"/>
   <fraction n="0.049" ref="aluminum"/>
   <fraction n="0.019" ref="iron"/>
   <fraction n="0.015" ref="potassium"/>
  </material>

  <material formula=" " name="ShotRock">
   <D value="@{[2.7*0.6]}" unit="g/cm3"/>
   <fraction n="0.438" ref="oxygen"/>
   <fraction n="0.257" ref="silicon"/>
   <fraction n="0.222" ref="sodium"/>
   <fraction n="0.049" ref="aluminum"/>
   <fraction n="0.019" ref="iron"/>
   <fraction n="0.015" ref="potassium"/>
  </material>

  <material formula=" " name="Dirt">
   <D value="1.7" unit="g/cm3"/>
   <fraction n="0.438" ref="oxygen"/>
   <fraction n="0.257" ref="silicon"/>
   <fraction n="0.222" ref="sodium"/>
   <fraction n="0.049" ref="aluminum"/>
   <fraction n="0.019" ref="iron"/>
   <fraction n="0.015" ref="potassium"/>
  </material>

  <material formula=" " name="Concrete">
   <D value="2.3" unit="g/cm3"/>
   <fraction n="0.530" ref="oxygen"/>
   <fraction n="0.335" ref="silicon"/>
   <fraction n="0.060" ref="calcium"/>
   <fraction n="0.015" ref="sodium"/>
   <fraction n="0.020" ref="iron"/>
   <fraction n="0.040" ref="aluminum"/>
  </material>

  <material formula="H2O" name="Water">
   <D value="1.0" unit="g/cm3"/>
   <fraction n="0.1119" ref="hydrogen"/>
   <fraction n="0.8881" ref="oxygen"/>
  </material>

  <material formula="Ti" name="Titanium">
   <D value="4.506" unit="g/cm3"/>
   <fraction n="1." ref="titanium"/>
  </material>

  <material name="TPB" formula="TPB">
   <D value="1.40" unit="g/cm3"/>
   <fraction n="1.0000" ref="argon"/>
  </material>

  <material name="Glass">
   <D value="2.74351" unit="g/cm3"/>
   <fraction n="0.600" ref="SiO2"/>
   <fraction n="0.118" ref="Al2O3"/>
   <fraction n="0.001" ref="Fe2O3"/>
   <fraction n="0.224" ref="CaO"/>
   <fraction n="0.034" ref="MgO"/>
   <fraction n="0.010" ref="Na2O"/>
   <fraction n="0.013" ref="TiO2"/>
  </material>

  <material name="Acrylic">
   <D value="1.19" unit="g/cm3"/>
   <fraction n="0.600" ref="carbon"/>
   <fraction n="0.320" ref="oxygen"/>
   <fraction n="0.080" ref="hydrogen"/>
  </material>

</materials>
EOF

close(MAT);
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++ gen_TPC ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


sub gen_TPC()
{

#constructs everything inside volTPC, namely
# (moving from left to right, or from +x to -x)
#  -volTPCPlaneU: with wires angled from vertical
#  -volTPCPlaneV: with wires angled from vertical 
#  -volTPCPlaneX: with vertical wires


# Create the TPC fragment file name,
# add file to list of output GDML fragments,
# and open it
    $TPC = "icarus_TPC" . $suffix . ".gdml";
    push (@gdmlFiles, $TPC);
    $TPC = ">" . $TPC;
    open(TPC) or die("Could not open file $TPC for writing");


# The standard XML prefix and starting the gdml
    print TPC <<EOF;
<?xml version='1.0'?>
<gdml>
EOF


# All the TPC solids save the wires.
print TPC <<EOF;
<solids>
    <box name="TPC" lunit="cm" 
      x="$TPC_x" 
      y="$TPC_y" 
      z="$TPC_z"/>
    <box name="TPCPlane" lunit="cm" 
      x="$TPCWirePlane_x" 
      y="$TPCWirePlane_y" 
      z="$TPCWirePlane_z"/>
    <box name="TPCActive" lunit="cm"
      x="$TPCActive_x"
      y="$TPCActive_y"
      z="$TPCActive_z"/>
EOF


#++++++++++++++++++++++++++++ Wire Solids ++++++++++++++++++++++++++++++

print TPC <<EOF;

    <tube name="TPCWireVert"
      rmax="@{[0.5*$TPCWireThickness]}"
      z="$TPCWirePlane_z"               
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

# Set number of wires to default to zero, when $wires_on = 0, for a low memory 
# version. But if $wires_on = 1, calculate the number of wires on each side of each
# plane to be used in the for loops

my $NumberCornerUWires = 0;
my $NumberSideUWires = 0;
my $NumberCommonUWires = 0;
my $NumberCornerVWires = 0;
my $NumberSideVWires = 0;
my $NumberCommonVWires = 0;
my $NumberVerticalWires = 0;

if ($wires_on == 1)
{
   # Number of wires in one corner
    $NumberCornerUWires = int( $TPCWirePlane_y/$UWire_ypitch );
    $NumberCornerVWires = int( $TPCWirePlane_y/$VWire_ypitch );

   # Total number of wires touching one vertical (longer) side
   # Note that the total number of wires per plane is this + another set of corner wires
    $NumberSideUWires = int( $TPCWirePlane_z/$UWire_zpitch );
    $NumberSideVWires = int( $TPCWirePlane_z/$VWire_zpitch );

   # Number of wires per side that aren't cut off by the corner
    $NumberCommonUWires = $NumberSideUWires - $NumberCornerUWires;
    $NumberCommonVWires = $NumberSideVWires - $NumberCornerVWires;

   # number of wires on the vertical plane
    $NumberVerticalWires = int( ($TPCWirePlane_y-$TPCWireThickness)/$XWirePitch );
}

# These XML comments throughout make the GDML file easier to navigate
print TPC <<EOF;

<!--+++++++++++++++++++ U Wire Solids ++++++++++++++++++++++-->

EOF

# The corner wires for the U plane
if ($wires_on==1) 
{
#VELE
    for ($i = 0; $i < $NumberCornerUWires; ++$i)
   {

	$length = ($i+1)*$UWire_zpitch / $SinUAngle;

   print TPC <<EOF;
    <tube name="TPCWireU$i"
      rmax="@{[0.5*$TPCWireThickness]}"
      z="@{[$length-0.06]}"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

    } #fine vele

   print TPC <<EOF;
    <tube name="TPCWireUCommon"
      rmax="@{[0.5*$TPCWireThickness]}"
      z="@{[$TPCWirePlane_y/$CosUAngle-0.06]}"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

} else { 

print TPC <<EOF;

        <!-- This GDML version has no wires and uses much less memory -->

EOF

}

print TPC <<EOF;


<!--+++++++++++++++++++ V Wire Solids ++++++++++++++++++++++-->


EOF

# The corner wires for the V plane
if ($wires_on==1) 
{
#VELE 
   for ($i = 0; $i < $NumberCornerVWires; ++$i)
    {
	$length = ($i+1)*$VWire_zpitch / $SinVAngle;

	print TPC <<EOF;

    <tube name="TPCWireV$i"
      rmax="@{[0.5*$TPCWireThickness]}"
      z="@{[$length-0.06]}"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>

EOF

   } #fine vele

    # The wire used many times in the middle of the V plane
    # Same subtraction as U common

   print TPC <<EOF;
    <tube name="TPCWireVCommon"
      rmax="@{[0.5*$TPCWireThickness]}"
      z="@{[$TPCWirePlane_y/$CosVAngle-0.06]}"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

} else { 

print TPC <<EOF;

                   <!-- no wires in this GDML -->

EOF

}


# Begin structure and create the vertical wire logical volume
print TPC <<EOF;
</solids>
<structure>
    <volume name="volTPCActive">
      <materialref ref="LAr"/>
      <solidref ref="TPCActive"/>
    </volume>


<!--+++++++++++++++++ Wire Logical Volumes ++++++++++++++++++++-->

EOF


if ($wires_on==1) 
{ 
  print TPC <<EOF;
    <volume name="volTPCWireVert">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCWireVert" />
    </volume>
EOF

  # Corner U wires logical volumes 
#VELE
  for ($i = 0; $i < $NumberCornerUWires; ++$i)
  {
  print TPC <<EOF;
    <volume name="volTPCWireU$i">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCWireU$i" />
    </volume>
EOF

  } #fine vele

  # Common U wire logical volume, referenced many times
  print TPC <<EOF;
    <volume name="volTPCWireUCommon">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCWireUCommon" />
    </volume>
EOF

  # Corner V wires logical volumes
#VELE
  for ($i = 0; $i < $NumberCornerVWires; ++$i)
  {
  print TPC <<EOF;
    <volume name="volTPCWireV$i">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCWireV$i" />
    </volume>
EOF

  } #fine vele

  # Common V wire logical volume, referenced many times
  print TPC <<EOF;
    <volume name="volTPCWireVCommon">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCWireVCommon" />
    </volume>
EOF

} else { 


print TPC <<EOF;

        <!-- This GDML version has no wires and uses much less memory -->

EOF

}


#+++++++++++++++++++++++++ Position physical wires ++++++++++++++++++++++++++

#            ++++++++++++++++++++++  U Plane  +++++++++++++++++++++++

# Create U plane logical volume
print TPC <<EOF;


<!--+++++++++++++++++++++ U Plane ++++++++++++++++++++++++-->


    <volume name="volTPCPlaneU">
      <materialref ref="LAr"/>
      <solidref ref="TPCPlane"/>
EOF

if ($wires_on==0)
{
print TPC <<EOF;

           <!-- no wires -->

EOF

} else {

#VELE
    for ($i = 0; $i < $NumberCornerUWires; ++$i)
    {

	if($i<$nSkipCornerUWires){ next; }
	my $ypos = - $TPCWirePlane_y/2 + (($i+1)*$UWire_ypitch)/2;
	my $zpos =   $TPCWirePlane_z/2 - (($i+1)*$UWire_zpitch)/2;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireU$i"/>
        <position name="posTPCWireU$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotationref ref="rPlusUAngleAboutX"/>
      </physvol>
EOF
    } #fine vele


    for ($i = 0; $i < $NumberSideUWires-$NumberCornerUWires; ++$i)
    {
	my $zpos =   $TPCWirePlane_z/2 
	           - ($NumberCornerUWires*$UWire_zpitch)/2
                   - ($i+1)*$UWire_zpitch;


print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireUCommon"/>
        <position name="posTPCWireU$i" unit="cm" x="0" y="0 " z="$zpos"/>
        <rotationref ref="rPlusUAngleAboutX"/>
      </physvol>
EOF
    }

#VELE
   for ($i = 0; $i < $NumberCornerUWires; ++$i)
   {
   
	if($i<$nSkipCornerUWires){ next; }
	my $ypos =   $TPCWirePlane_y/2 - (($i+1)*$UWire_ypitch)/2;
	my $zpos = - $TPCWirePlane_z/2 + (($i+1)*$UWire_zpitch)/2;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireU$i"/>
        <position name="posTPCWireU$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotationref ref="rPlusUAngleAboutX"/>
      </physvol>
EOF
    }#fine vele

} #ends else


#            ++++++++++++++++++++++  V Plane  +++++++++++++++++++++++

# End U plane and create V plane logical volume
print TPC <<EOF;
    </volume>




<!--+++++++++++++++++++++ V Plane ++++++++++++++++++++++++-->


    <volume name="volTPCPlaneV">
      <materialref ref="LAr"/>
      <solidref ref="TPCPlane"/>
EOF

if ($wires_on==0)
{
print TPC <<EOF;

           <!-- no wires -->

EOF

} else {

#VELE
   for ($i = 0; $i < $NumberCornerVWires; ++$i)
    {

	if($i<$nSkipCornerVWires){ next; }
	my $ypos =   $TPCWirePlane_y/2 - (($i+1)*$VWire_ypitch)/2;
	my $zpos =   $TPCWirePlane_z/2 - (($i+1)*$VWire_zpitch)/2;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireV$i"/>
        <position name="posTPCWireV$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotationref ref="rPlusVAngleAboutX"/>
      </physvol>
EOF
    }

    for ($i = 0; $i < $NumberSideVWires-$NumberCornerVWires; ++$i)
    {
	my $zpos =   $TPCWirePlane_z/2 
	           - ($NumberCornerVWires*$VWire_zpitch)/2
                   - ($i+1)*$VWire_zpitch;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireVCommon"/>
        <position name="posTPCWireV$i" unit="cm" x="0" y="0 " z="$zpos"/>
        <rotationref ref="rPlusVAngleAboutX"/>
      </physvol>
EOF
    }

#VELE
    for ($i = 0; $i < $NumberCornerVWires; ++$i)
    {
   
	if($i<$nSkipCornerVWires){ next; }
	my $ypos = - $TPCWirePlane_y/2 + (($i+1)*$VWire_ypitch)/2;
	my $zpos = - $TPCWirePlane_z/2 + (($i+1)*$VWire_zpitch)/2;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireV$i"/>
        <position name="posTPCWireV$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotationref ref="rPlusVAngleAboutX"/>
     </physvol>
EOF
    }#fine vele


} #ends else



#            ++++++++++++++++++++++  X Plane  +++++++++++++++++++++++

# End V plane and create X plane logical volume
print TPC <<EOF;
    </volume>



<!--+++++++++++++++++++++ X Plane ++++++++++++++++++++++++-->


    <volume name="volTPCPlaneX">
      <materialref ref="LAr"/>
      <solidref ref="TPCPlane"/>
EOF

if ($wires_on==0)
{
print TPC <<EOF;

           <!-- no wires -->

EOF

} else {


for ($i=0; $i<$NumberVerticalWires; ++$i)
{

my $ypos = (-0.5*$TPCWirePlane_y)+$TPCWireThickness/2+$XWirePitch*($i+0.5);

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireVert"/>
        <position name="posTPCWireX$i" unit="cm" x="0" y="$ypos " z="0"/>
        <rotationref ref="rIdentity" />
      </physvol>
EOF

}

} #ends else

print TPC <<EOF;
    </volume>
EOF

#+++++++++++++++++++++ Position physical wires Above +++++++++++++++++++++

#originale
# my $VolU_x = (-$TPC_x/2)+3*$WirePlaneSpacing;
# my $VolV_x = (-$TPC_x/2)+2*$WirePlaneSpacing;
# my $VolX_x = (-$TPC_x/2)+$WirePlaneSpacing;

my $VolU_x = (-$TPC_x/2)+3*$WirePlaneSpacing;
my $VolV_x = (-$TPC_x/2)+2*$WirePlaneSpacing;
my $VolX_x = (-$TPC_x/2)+$WirePlaneSpacing;



#wrap up the TPC file
print TPC <<EOF;
    <volume name="volTPC">
      <materialref ref="LAr" />
      <solidref ref="TPC" />
     <physvol>
       <volumeref ref="volTPCPlaneU" />
       <position name="posTPCPlaneU" unit="cm" x="$VolU_x" y="0" z="0" />
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneV" />
       <position name="posTPCPlaneV" unit="cm" x="$VolV_x" y="0" z="0" />
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneX" />
       <position name="posTPCPlaneX" unit="cm" x="$VolX_x" y="0" z="0" />
     </physvol>
     <physvol>
       <volumeref ref="volTPCActive"/>
       <positionref ref="posActiveInTPC"/>
     </physvol>
    </volume>
</structure>
</gdml>
EOF

    close(GDML);

} #end of sub gen_TPC


##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
sub gen_pmt {

    $PMT = "icarus_pmt" . $suffix . ".gdml";
    push (@gdmlFiles, $PMT); # Add file to list of GDML fragments
    $PMT = ">" . $PMT;
    open(PMT) or die("Could not open file $PMT for writing");

# The standard XML prefix and starting the gdml
    print PMT <<EOF;
<?xml version='1.0'?>
<gdml>
EOF

	print PMT <<EOF;
<solids>
 <tube name="PMTVolume"
  rmax="10.16"
  z="2.54"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
</solids>
<structure>
 <volume name="volOpDetSensitive">
  <materialref ref="Glass"/>
  <solidref ref="PMTVolume"/>
 </volume>
</structure>
EOF

#Close standard XML file
print PMT <<EOF;
</gdml>
EOF

}
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_Cryostat +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Cryostat()
{

# Create the cryostat fragment file name,
# add file to list of output GDML fragments,
# and open it
    $CRYO = "icarus_Cryostat" . $suffix . ".gdml";
    push (@gdmlFiles, $CRYO);
    $CRYO = ">" . $CRYO;
    open(CRYO) or die("Could not open file $CRYO for writing");


# The standard XML prefix and starting the gdml
    print CRYO <<EOF;
<?xml version='1.0'?>
<gdml>
EOF


# All the cryostat solids.
print CRYO <<EOF;
<solids>

    <box name="Cryostat" lunit="cm" 
      x="$Cryostat_x" 
      y="$Cryostat_y" 
      z="$Cryostat_z"/>
    <box name="ArgonInterior" lunit="cm" 
      x="$LAr_x"
      y="$LAr_y"
      z="$LAr_z"/>
    <box name="GaseousArgon" lunit="cm" 
      x="$LAr_x"
      y="$GaseousAr_y"
      z="$LAr_z"/>
    <subtraction name="SteelShell">
      <first ref="Cryostat"/>
      <second ref="ArgonInterior"/>
    </subtraction>

    <box name="Cathode" lunit="cm"
      x="$CPA_x"
      y="$TPC_y"
      z="$TPC_z"/>

</solids>
EOF

# Cryostat structure
print CRYO <<EOF;
<structure>
    <volume name="volSteelShell">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="SteelShell" />
    </volume>
    <volume name="volGaseousArgon">
      <materialref ref="ArGas"/>
      <solidref ref="GaseousArgon"/>
    </volume>

    <volume name="volCathode">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="Cathode" />
    </volume>


    <volume name="volCryostat">
      <materialref ref="LAr" />
      <solidref ref="Cryostat" />
      <physvol>
        <volumeref ref="volGaseousArgon"/>
        <position name="posGaseousArgon" unit="cm" x="0" y="@{[$LAr_y/2-$GaseousAr_y/2]}" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volSteelShell"/>
        <position name="posSteelShell" unit="cm" x="0" y="0" z="0"/>
      </physvol>

      <physvol>
        <volumeref ref="volTPC"/>
        <positionref ref="posTPC0inCryo"/>
	<rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volCathode" />
        <positionref ref="posCathode"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPC"/>
        <positionref ref="posTPC1inCryo"/>
	<rotationref ref="rPlus180AboutY"/>
      </physvol>

EOF
############################################################################################
############################################################################################
############################################################################################
############################################################################################
### @pmt_pos = ( ' x="-74.8"  y="0"  z="0"',
###             ' x=" 74.8" y="0" z="0"');

#$PMT_x = (-$TPC_x/2)+3*$WirePlaneSpacing;

$PMT_x0 = - $TPC_x - 10; 
@pmt_pos0 = read_pmt_pos("disposizione_caso.txt", $PMT_x0);
$Num_PMTs0 = @pmt_pos0;

$PMT_x1 =   $TPC_x + 10;
@pmt_pos1 = read_pmt_pos("disposizione_caso.txt", $PMT_x1);
$Num_PMTs1 = @pmt_pos1;


    for ( $i=0; $i<$Num_PMTs0; ++$i ){
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volOpDetSensitive"/>
   <position name="posPMT0$i" unit="cm" @pmt_pos0[$i]/>
   <rotationref ref="rPlus90AboutY"/>
  </physvol>
EOF
    }

    for ( $i=0; $i<$Num_PMTs1; ++$i ){
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volOpDetSensitive"/>
   <position name="posPMT1$i" unit="cm" @pmt_pos1[$i]/>
   <rotationref ref="rPlus90AboutY"/>
  </physvol>
EOF
    }

#    for ( $i=0; $i<$Num_PMTs; ++$i ){
#      print CRYO <<EOF;
#  <physvol>
#   <volumeref ref="volPMTBig"/>
#   <position name="posPMT$i" unit="cm" @pmt_pos[$i]/>
#   <rotationref ref="rPlus90AboutY"/>
#  </physvol>
#EOF
#    }

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

print CRYO <<EOF;
    </volume>
</structure>
</gdml>
EOF

close(CRYO);
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++ gen_Enclosure +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Enclosure()
{

# Create the detector enclosure fragment file name,
# add file to list of output GDML fragments,
# and open it
    $ENCL = "icarus_DetEnclosure" . $suffix . ".gdml";
    push (@gdmlFiles, $ENCL);
    $ENCL = ">" . $ENCL;
    open(ENCL) or die("Could not open file $ENCL for writing");


# The standard XML prefix and starting the gdml
    print ENCL <<EOF;
<?xml version='1.0'?>
<gdml>
EOF


# All the detector enclosure solids.
print ENCL <<EOF;
<solids>

    <box name="DetEnclosure" lunit="cm" 
      x="$DetEnc_x"
      y="$DetEnc_y"
      z="$DetEnc_z"/>

</solids>
EOF



# Detector enclosure structure
    print ENCL <<EOF;
<structure>

    <volume name="volDetEnclosure">
      <materialref ref="Concrete"/>
      <solidref ref="DetEnclosure"/>

      <physvol>
        <volumeref ref="volCryostat"/>
        <positionref ref="posCryo1InDetEnc"/>
      </physvol>

      <physvol>
        <volumeref ref="volCryostat"/>
        <positionref ref="posCryo2InDetEnc"/>
      </physvol>
    </volume>
EOF

print ENCL <<EOF;
</structure>
</gdml>
EOF

close(ENCL);
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++ gen_World +++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_World()
{

# Create the WORLD fragment file name,
# add file to list of output GDML fragments,
# and open it
    $WORLD = "icarus_World" . $suffix . ".gdml";
    push (@gdmlFiles, $WORLD);
    $WORLD = ">" . $WORLD;
    open(WORLD) or die("Could not open file $WORLD for writing");


# The standard XML prefix and starting the gdml
    print WORLD <<EOF;
<?xml version='1.0'?>
<gdml>
EOF


# All the World solids.
print WORLD <<EOF;
<solids>
    <box name="World" lunit="cm" 
      x="$World_x" 
      y="$World_y" 
      z="$World_z"/>
</solids>
EOF

# World structure
print WORLD <<EOF;
<structure>
    <volume name="volWorld" >
      <materialref ref="Air"/>
      <solidref ref="World"/>

      <physvol>
        <volumeref ref="volDetEnclosure"/>
	<positionref ref="posDetEncInWorld"/>
      </physvol>

    </volume>
</structure>
</gdml>
EOF

# make_gdml.pl will take care of <setup/>

close(WORLD);
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++ write_fragments ++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub write_fragments()
{
   # This subroutine creates an XML file that summarizes the the subfiles output
   # by the other sub routines - it is the input file for make_gdml.pl which will
   # give the final desired GDML file. Specify its name with the output option.
   # (you can change the name when running make_gdml)

   # This code is taken straigh from the similar MicroBooNE generate script, Thank you.

    if ( ! defined $output )
    {
	$output = "-"; # write to STDOUT 
    }

    # Set up the output file.
    $OUTPUT = ">" . $output;
    open(OUTPUT) or die("Could not open file $OUTPUT");

    print OUTPUT <<EOF;
<?xml version='1.0'?>

<!-- Input to Geometry/gdml/make_gdml.pl; define the GDML fragments
     that will be zipped together to create a detector description. 
     -->

<config>

   <constantfiles>

      <!-- These files contain GDML <constant></constant>
           blocks. They are read in separately, so they can be
           interpreted into the remaining GDML. See make_gdml.pl for
           more information. 
	   -->
	   
EOF

    foreach $filename (@defFiles)
    {
	print OUTPUT <<EOF;
      <filename> $filename </filename>
EOF
    }

    print OUTPUT <<EOF;

   </constantfiles>

   <gdmlfiles>

      <!-- The GDML file fragments to be zipped together. -->

EOF

    foreach $filename (@gdmlFiles)
    {
	print OUTPUT <<EOF;
      <filename> $filename </filename>
EOF
    }

    print OUTPUT <<EOF;

   </gdmlfiles>

</config>
EOF

    close(OUTPUT);
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub read_pmt_pos {

  $pmt_x = $_[1];

  $PMT_pos_file = $_[0];
  open(PMTPOS, $PMT_pos_file) or die("Could not open file $PMT_pos_file.");

  @pmt_pos = ();

  foreach $line (<PMTPOS>) {

    @coord = split(/\s/, $line);

    $string = " x=\" $pmt_x\" y=\"$coord[1]\" z=\"$coord[0]\"";
    push(@pmt_pos, $string);

  }

  close(PMTPOS);

  return @pmt_pos;
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#copy from generate_gdml_test.pl with different parts of PMTs

sub gen_pmtBig {

    $PMTBig = "micro-pmtdef" . $suffix . ".gdml";
    push (@gdmlFiles, $PMTBig); # Add file to list of GDML fragments
    $PMTBig = ">" . $PMTBig;
    open(PMTBig) or die("Could not open file $PMTBig for writing");

 print PMTBig <<EOF;
<?xml version='1.0'?>
<gdml>
EOF

	print PMTBig <<EOF;
<solids>
 <tube name="PMTVolumeBig"
  rmax="15"
  z="27.50"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMT_TPBCoating"
  rmax="15"
  z="0.01"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMT_AcrylicPlate"
  rmax="15"
  z="0.2"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMT_Stalk"
  rmax="3.75"
  z="7.5"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMT_SteelBase"
  rmax="15"
  z="3.75"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMT_Underside"
  rmax="10"
  z="6.25"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMT_Lens"
  rmax="10"
  z="6.25"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
</solids>
<structure>
 <volume name="vol_PMT_TPBCoating">
  <materialref ref="TPB"/>
  <solidref ref="PMT_TPBCoating"/>
 </volume>
 <volume name="vol_PMT_AcrylicPlate">
  <materialref ref="Acrylic"/>
  <solidref ref="PMT_AcrylicPlate"/>
 </volume>
 <volume name="vol_PMT_Stalk">
  <materialref ref="Glass"/>
  <solidref ref="PMT_Stalk"/>
 </volume>
 <volume name="vol_PMT_SteelBase">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="PMT_SteelBase"/>
 </volume>
 <volume name="vol_PMT_Underside">
  <materialref ref="Glass"/>
  <solidref ref="PMT_Underside"/>
 </volume>
 <volume name="volOpDetSensitive">
  <materialref ref="LAr"/>
  <solidref ref="PMT_Lens"/>
 </volume>
 <volume name="volPMTBig">
  <materialref ref="LAr"/>
  <solidref ref="PMTVolumeBig"/>
  <physvol>
   <volumeref ref="vol_PMT_TPBCoating"/>
   <position name="pos_PMT_TPBCoating" unit="cm" x="0" y="0" z="13.97"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_AcrylicPlate"/>
   <position name="pos_PMT_AcrylicPlate" unit="cm" x="0" y="0" z="13.86"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_Stalk"/>
   <position name="pos_PMT_Stalk" unit="cm" x="0" y="0" z="-6.25"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_SteelBase"/>
   <position name="pos_PMT_SteelBase" unit="cm" x="0" y="0" z="-12"/>
  </physvol>
  <physvol>
   <volumeref ref="volOpDetSensitive"/>
   <position name="pos_PMT_Lens" unit="cm" x="0" y="0" z="-3.81"/>*  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_Underside"/>
   <position name="pos_PMT_Underside" unit="cm" x="0" y="0" z="-3.81"/>
  </physvol>
 </volume>
</structure>
EOF

#Close standard XML file
print PMTBig <<EOF;
</gdml>
EOF

}








