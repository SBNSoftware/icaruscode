#!/usr/bin/perl


# Each subroutine generates a fragment GDML file, and the last subroutine
# creates an XML file that make_gdml.pl will use to appropriately arrange
# the fragment GDML files to create the final desired ICARUS GDML file, 
# to be named by make_gdml output command

# If you are playing with different geometries, you can use the
# suffix command to help organize your work.

use Math::Trig;
use Getopt::Long;
use Math::BigFloat;
Math::BigFloat->precision(-16);
#Math::BigFloat->accuracy(16);

GetOptions( "help|h" => \$help,
	    "suffix|s:s" => \$suffix,
	    "output|o:s" => \$output,
	    "concrete|c:s" => \$thickness_over,
	    "wires|w:s" => \$wires,
	    "vetocrt|v:s" => \$crt);

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

# set wires on to be the default, unless given an input by the user:  1=on, 0=off
if (defined $wires)
{
    $wires_on = $wires;

}

else { $wires_on = 1;}   # 1=on, 0=off


#set crt on or off: 1=on, 0=off
if (defined $crt)
{
    $crt_on = $crt;
}

else { $crt_on = 1;} # 1=on, 0=off

#set thickness of the concrete overburden. Remmeber: the dafaul value is 300 cm and 0 means no everburden
if (defined $thickness_over)
{
	$concrete_on = $thickness_over; 
}

else { $concrete_on = 300;}


#-------Definitions of all variables: unit= cm

$inch = 2.54;

##################################################################
##################### wire plane parameters ######################

$YWirePitch             =   0.3;
$UWirePitch             =   0.3;
$VWirePitch             =   0.3;

#respect to the z axes
#$YAngle                 =   0;
$UAngle                 =   60;
$VAngle			=   60;

$SinUAngle              =   sin( deg2rad($UAngle) );
$CosUAngle              =   cos( deg2rad($UAngle) );
$TanUAngle              =   tan( deg2rad($UAngle) );

$SinVAngle              =   sin( deg2rad($VAngle) );
$CosVAngle              =   cos( deg2rad($VAngle) );
$TanVAngle              =   tan( deg2rad($VAngle) );

#$UWireCornerInt_y       =   $UWirePitch * $CosUAngle;
#$UWireCornerInt_z       =   $UWirePitch * $SinUAngle;
$UWire_ypitch           =   $UWirePitch / $CosUAngle;
$UWire_zpitch           =   $UWirePitch / $SinUAngle;

#$VWireCornerInt_y       =   $VWirePitch * $CosVAngle;
#$VWireCornerInt_z       =   $VWirePitch * $SinVAngle;
$VWire_ypitch           =   $VWirePitch / $CosVAngle;
$VWire_zpitch           =   $VWirePitch / $SinVAngle;


$TPCWireThickness       =   0.015; #wire diameter
$TPCWireRadius = $TPCWireThickness /2;  #wire radius

$CommonWireLength       =   (1056 * $YWirePitch) / $SinUAngle; #364.90;
#print("CommonWireLength: $CommonWireLength \n");

###########################################################################
########################### spacing parameters ############################


$CPA_x                 =     0.15;  #cathode plane 1.5 mm
$WirePlaneSpacing      =     0.3;   # center to center
$MaxDrift              =     148.2; #drift length in LAr at cryogenic temperature
 
#Cryostat space with LAr outside of entire fiducial volume
$SpaceWirePlToWall     =     31.8; 
$SpaceWirePlToWirePl   =     85; # edge to edge, like all of these (it was in the original perl script)
#$SpaceTPCToFloor       =     37; # from the article
#$SpaceTPCToTopLAr      =     30.5;  
$SpaceTPCToFloor       =     36; # from the article
$SpaceTPCToTopLAr      =     29.5; 
$UpstreamLArPadding    =     82.50; # from the article
$DownstreamLArPadding  =     82.50; # from the article

##########################################################
######## Other Cathode relevant parameters ###############

$CPA_epsi   = 0.001; #to avoid overlaps
$CPA_y      = 320; #Cathode height.
$CPA_z      = 1800; #Cathode length.
$CPAStrip_x = $CPA_x-$CPA_epsi; #Cathode strip thickness.
$CPAStrip_y = 2.1; #Cathode strip height.
$CPAStrip_pitch = 5; #Cathode strip pitch.

##############################################################
############## Cryo and TPC relevant dimensions  #############

#Active LAr volume
$TPCActive_x    =     $MaxDrift;
$TPCActive_y    =     $CommonWireLength * $SinUAngle + 0.02; #316.0; 
$TPCActive_z    =     1795.5;

#print("TPCActive_x: $TPCActive_x, TPCActive_y: $TPCActive_y, TPCActive_z: $TPCActive_z ");

    #Reset the active length and the common wire length
$TPCActive_z    =     $CommonWireLength * $CosUAngle + (4640 - 1) * $UWire_zpitch;

#print(" - change with wire pitch: $UWire_zpitch to $TPCActive_z \n");

# TPCWirePlane dimensions
$TPCWirePlane_x     =       2*$TPCWireThickness; #ATTENTION: tentative to not overlap planes: is this correct?Should be the right thickness of the planes!
$TPCWirePlane_y     =       $TPCActive_y; 
$TPCWirePlane_z     =       $TPCActive_z;

#print("TPCWirePlane_x: $TPCWirePlane_x, TPCWirePlane_y: $TPCWirePlane_y, TPCWirePlane_z: $TPCWirePlane_z \n");

##################################################################
#Dimension of the TPC (active+passive volume)
$TPC_x    =     $MaxDrift+ 6*$TPCWireThickness + 3*$WirePlaneSpacing + $CPA_x;
#$TPC_x    =     150.0;

$DeltaLUCorner = $UWirePitch/($SinUAngle*$CosUAngle); #this is the Delta L for the corner wire length
$DeltaLVCorner = $VWirePitch/($SinVAngle*$CosVAngle);

#print(" ******** DELTA CORENR ***** \n");
#print(" $DeltaLUCorner  $DeltaLVCorner \n");

$LAr_x    =     $CPA_x 
              + 2*($TPC_x + $SpaceWirePlToWall);
$LAr_y    =     $TPCActive_y 
              + $SpaceTPCToFloor 
              + $SpaceTPCToTopLAr;
$LAr_z    =     $TPCActive_z
              + $UpstreamLArPadding 
              + $DownstreamLArPadding;

#$TPC_y    =     388.0; #From official drawings
#$TPC_z    =     1969.9; #From official drawings
$TPC_y    =     $LAr_y;
$TPC_z    =     $LAr_z;
    
$AlumThickness		=	19;    
$GaseousAr_y            =       6.5;
$CryoDist 		=	20;

$Cryostat_x = $LAr_x + 2*$AlumThickness ; 
#$Cryostat_y = $LAr_y + 2*$AlumThickness + $GaseousAr_y ;
#$Cryostat_z = $LAr_z + 2*$AlumThickness ;
#$Cryostat_y = $TPC_y + 2*$AlumThickness ;
$Cryostat_y = $TPC_y + $GaseousAr_y + 2*$AlumThickness ;
$Cryostat_z = $TPC_z + 2*$AlumThickness ;

#$LAr_x_orig =   4*($TPC_x) + 2*($CPA_x + $SpaceWirePlToWall) + $SpaceWirePlToWirePl;

$LAr_x_orig = 2*$LAr_x + 2*$AlumThickness + $CryoDist ;	    #for total positioning

$Cryostat_x_orig = $LAr_x_orig + 2*$AlumThickness ;               #for total positioning


$TPCinCryo_x[0]     =      - $TPC_x/2 - $CPA_x/2;
#$TPCinCryo_x[1]     =      - $LAr_x/2 + $SpaceWirePlToWall + 1.5*($TPC_x) + $CPA_x ;
#$TPCinCryo_x[2]     =        $LAr_x/2 - $SpaceWirePlToWall - 1.5*($TPC_x) - $CPA_x ;
$TPCinCryo_x[1]     =        $TPC_x/2 + $CPA_x/2;

$posCat_x      =      0;

#$posRightCat_x      =      - $Cryostat_x/2 + 1.0*($TPC_x) + $SpaceWirePlToWall + $CPA_x/2 ;
#$posLeftCat_x       =        $Cryostat_x/2 - 1.0*($TPC_x) - $SpaceWirePlToWall - $CPA_x/2 ;

#$TPCinCryo_y        =      - $Cryostat_y/2 + $TPC_y/2 + $SpaceTPCToFloor;  
#$TPCinCryo_z        =      - $Cryostat_z/2 + $TPC_z/2 + $UpstreamLArPadding;  
#$TPCinCryo_y        =      0.;    
$TPCinCryo_y        =      -$GaseousAr_y/2;
$TPCinCryo_z        =      0.; 

$posTPCActive_x     = 0.570000000000007;
$posTPCActive_y     = 0;
$posTPCActive_z     = 0;

##################################################################
#am ##########  Steel Mechanical Structure Parameters ############
##################################################################

  $struct_width = 10.;     # width of the structure box segment;
  $ext_struct_x = $struct_width;   # external width of structure segment;
  $ext_struct_y = 341.8;   # external heigh of structure segment;
  $ext_struct_z = 199.532; # external length of structure segment;
  $int_struct_x = $ext_struct_x;  # internal width of structure segment;
  $int_struct_y = $ext_struct_y - 2*$struct_width;   # internal heigh of structure segment;
  $int_struct_z = $ext_struct_z - 2*$struct_width;   # internal length of structure segment;

##################################################################
############## PMTs relevant parameters  #########################

$NumberPMT = 90;
$PMTthickness = 0.5;
$PMTradius = 4*$inch;
$PMTradiusOuter = 4*$inch; #10.16 cm
$PMTradiusInner = $PMTradiusOuter - $PMTthickness;
$PMTGrid_offset = 2; # from GLR communication
$PMTx = 0.557*$PMTradiusOuter; # x semi-axis of PMT
#$PMTtube_x = 10.3; #length of the glass tube
$PMTtubeRmax = 4.225; #outer radius of PMT glass tube;
$PMTtubeRmin = 3.725; #inner radius of PMT glass tube;$PMTtube_x = 2*$PMTtubeRmax;
$PMTtube_x = 2*$PMTtubeRmax ;
#$PMTxCut = ($PMTx*$PMTx - $PMTtubeRmax*$PMTtubeRmax)^0.5;
$PMTPlane_x = 2*$PMTradius+$PMTGrid_offset+$PMTtube_x;
#$PMTPlane_y = $TPCActive_y ; 
$PMTPlane_y = $ext_struct_y + 6.;
#$PMTPlane_z = $TPCActive_z ; 
$PMTPlane_z = 9*$ext_struct_z + 2*struct_width + 6.;
$PMTWiresOffset = 0.5 ; # space between PMTs and wire planes 

##################################################################
######################### CRT parameters  ########################

#defined in seperate python script "gen_crt_frags.py"
#N.B. warm vessel origin, size currently hardcoded - changes
# in these values here need to be copied over (will fix eventually..)

##################################################################
############## DetEnc and World relevant parameters  #############

#Foam thickness and Warm vessel thickness
$FoamPadding            =       60; 
$WarmVesselThickness    =       27.4;
$WarmVesselThickness_x  =       21.8;
$SpaceTop               =       22.4; #Space between cryostat and top thermal insulation
$SpaceUpstream          =       65.55; #Space between cryostat and upstream thermal insulation
$SpaceDownstream        =       15.55; #Space between cryostat and downstream thermal insulation
$SpaceBottom            =       18.0; #Space between cryostat and bottom thermal insulation 

#Thermal insulation surrounding the two cryostats
$ThermIns_x             =       2*$Cryostat_x+2*$FoamPadding + 3*$CryoDist;
$ThermIns_y	        =	$Cryostat_y+2*$FoamPadding+$SpaceTop+$SpaceBottom; 
$ThermIns_z             =       $Cryostat_z+2*$FoamPadding+$SpaceUpstream+$SpaceDownstream;

#$InDetEncl_y            =       -154.995; #0; #-$WarmVessel_y/6; #-($WarmVessel_y + $CRT_tot_y)/6; original value
#$InDetEncl_y            =       -154.995; #0; #-$WarmVessel_y/6; #-($WarmVessel_y + $CRT_tot_y)/6; +15cm

$ThermInsInDetEncl_x    =       0;
#$ThermInsInDetEncl_y    =       $InDetEncl_y;
$ThermInsInDetEncl_z    =	0;

#Warm Vessel surrounding the thermal insulation
#$WarmVessel_x	        =	$ThermIns_x+2*$WarmVesselThickness;
$WarmVessel_x	        =	$ThermIns_x+2*$WarmVesselThickness_x;
$WarmVessel_y	        =	$ThermIns_y+2*$WarmVesselThickness; 
$WarmVessel_z           =       $ThermIns_z+2*$WarmVesselThickness;

$WarmVesselInDetEncl_x  =	0;
#$WarmVesselInDetEncl_y  =       $InDetEncl_y;
$WarmVesselInDetEncl_z  =	0;

#$MoreDetEnc_y = 60.0; # More space in height to cope with few extrusions
#Big detector Enclosure to contain detector + CRT.
#$DetEnc_x = 1235.96;#$WarmVessel_x + $CRT_tot_x; 
#$DetEnc_y = 963.37 + $MoreDetEnc_y;#$WarmVessel_y + $CRT_tot_y; 
#$DetEnc_z = 2709.56;#$WarmVessel_z + $CRT_tot_z; 


#Cryostat respect to the warm vessel
$Cryo1InWarmVessel_x     =     -$Cryostat_x/2 - $CryoDist ;  #-$Cryostat_x/2 - $CryoDist/2 ;
$Cryo2InWarmVessel_x     =      $Cryostat_x/2 + $CryoDist;  #$Cryostat_x/2 + $CryoDist/2 ;
#$CryoInWarmVessel_y     =       -$WarmVessel_y/2 + $ConcretePadding + $Cryostat_y/2; #in original way
#$CryoInWarmVessel_y     =       -$WarmVessel_y/2 + $FoamPadding + $WarmVesselThickness + $Cryostat_y/2 -154.995; #(-$WarmVessel_y/2 + $TotalPadding + $Cryostat_y/2)- ($WarmVessel_y)/6 ;
$CryoInWarmVessel_z     =       0;

#Original Origin point found in the Larsoft perl script
  # We want the world origin to be at the very front of the fiducial volume.
  # move it to the front of the enclosure, then back it up through the concrete/foam, 
  # then through the Cryostat shell, then through the upstream dead LAr (including the
  # dead LAr on the edge of the TPC, but this is covered in $UpstreamLArPadding).
  # This is to be added to the z position of every volume in volWorld

#$OriginZSet             =       $WarmVessel_z/2  - $TotalPadding  - $SteelThickness  - $UpstreamLArPadding;
#$OriginYSet             =       $DetEnWarmVessel_y/2 - $TotalPadding - $SteelThickness - $SpaceTPCToFloor - $TPC_y/2;
##$OriginXSet             =       $LAr_x_orig/2 - $SpaceWirePlToWall  - 3*$WirePlaneSpacing - $TPCWirePlane_x;
#$OriginXSet             =       $WarmVessel_x/2- $TotalPadding - $SteelThickness;

#Marta: I would like that the origin is in the middle of the two cryostats. It corresponds to the world center. It do not correspond to the experimental hall centre.
$OriginXSet = 0;
$OriginYSet = 0; 
$OriginZSet = 0;


#Experimental hall
$Building_y = 1040.0 ; #part of the building outside
$ExpHall_y = 1170.0 + 80;  #building underground

$Hall_x = 1890.0;
$Hall_y = $Building_y + $ExpHall_y ;
$Hall_z = 3870.0;

$HallWallThicnekss = 60;
$Overburden_yDefault = 300; #for positioning of the dirt volume: it not depends on the concrete dimensions

#True overburden dimensions
#$Overburden_x	    =	1830;  
#$Overburden_y	    =   $concrete_on;
#$Overburden_z 	    =	2930;

#Overburden dimensions to include also the floor
$Overburden_x	    =	1830.0 ;  #Hall_x - thickness of the wall
$Overburden_y	    =   $concrete_on;
$Overburden_z 	    =	3810.0; #Hall_z - thickness of the wall


#Overburden
$posOverburden_x 	= 	$OriginXSet ;
#$posOverburden_y 	=       $WarmVessel_y/2 + $Overburden_yDefault/2 + (820 - $WarmVessel_y/2) ; #warm vessel_y + overburden_y + distance between overburden and warm vessel 820 (taken from the drawings)
$posOverburden_z 	=	$OriginZSet;

#DetectorEnclosure (centered between expHall walls, between expHall floor and bottom of overburden)
$ExpHall_VertSpace = $ExpHall_y - $HallWallThicnekss - $Overburden_y/3;
$DetEnc_pad = 0.1;
$DetEnc_x = 1245;
$DetEnc_y = $ExpHall_VertSpace - 2*$DetEnc_pad;
$DetEnc_z = 3200;
$WarmVessel_FootHeight = 6.0; # heigth of support feet from bottom of support structure
$WarmVessel_FloorSpace = 10.16; # space between WV feet and concrete floor due to concrete islands, grout, padding
$WarmVessel_CenterToFloor = $WarmVessel_FootHeight + $WarmVessel_FloorSpace + $WarmVessel_y/2;
#$DetEncl_yOffset = $ExpHall_VertSpace/2 - $WarmVessel_CenterToFloor;
#$ThermInsInDetEncl_y = -1*$DetEncl_yOffset;
#$WarmVesselInDetEncl_y = $ThermInsInDetEncl_y;
#World
$World_x            =       1e4;	#Originally was 2*$DetEnc_x;
$World_y            =       1e4;	#Originally was 2*$DetEnc_y;
$World_z            =       1e4;	#Originally was 2*$DetEnc_z;


#Ground Level
#$Ground_y = -$World_y/4+$ExpHall_y/2-$Overburden_y/6 + 0.5*($World_y/2 + $ExpHall_y - $Overburden_y/3);
#$Ground_y = $ExpHall_y - $Overburden_y/3 ;
$Ground_y = 780.0;# - 40.0; #from detector building drawing: distance between beamline and ground level
$DetEncl_yOffset = $Ground_y - $Overburden_y/3 - $ExpHall_VertSpace/2 -$DetEnc_pad;
$ThermInsInDetEncl_y = -1*$DetEncl_yOffset;
$WarmVesselInDetEncl_y = $ThermInsInDetEncl_y;
$CryoInWarmVessel_y = $WarmVesselInDetEncl_y -$WarmVessel_y/2 + $FoamPadding + $WarmVesselThickness + $Cryostat_y/2;

$CRTSHELL_WV_OFFSET_Y = 15;
$CRTSHELL_WV_OFFSET_Z = 164.5865; #v08_22_00 had 118.2315;
#$CryoInWarmVessel_y = $WarmVesselInDetEncl_y -$WarmVessel_y/2 + $FoamPadding + $WarmVesselThickness + $Cryostat_y/2 + $GaseousAr_y/2;
#+++++++++++++++++++++++++ End defining variables ++++++++++++++++++++++++++

# run the sub routines that generate the fragments

gen_Define(); 	 # generates definitions at beginning of GDML
gen_Materials(); # generates materials to be used

gen_RaceTracks(); # generates the geometry of race-tracks

gen_TPC();	 # generates wires, wire planes, and puts them in volTPC
	         # This is the bulk of the code, and has zero wires option


gen_PMT();	 #generates PMTs

gen_Mech_Structure(); # generates the geometry of the TPC mechanical structure

gen_Cryostat();	 # places  volTPC,
		 # half rotated 180 about Y
if ($crt_on==1) {gen_CRT();}	 # places CRT: CERN modules top, eves; MINOS modules sides, DC bottom
gen_Enclosure(); # places two cryostats and warm vessel

gen_World();	 # places the enclosure in the experimental hall


write_fragments(); # writes the XML input for make_gdml.pl
			# which zips together the final GDML
exit;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++ usage +++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub usage()
{
    print "Usage: $0 [-h|--help] [-o|--output <fragments-file>] [-s|--suffix <string>]";
    print " [-c|--concrete <double>] [-w|--wires <wire or no wire geometry>] [-v <crt or no crt geometry>] \n";
    print "       if -o is omitted, output goes to STDOUT; <fragments-file.xml> is input to make_gdml.pl\n";
    print "       -s <string> appends the string to the file names; useful for multiple detector versions\n";
    print "       -c <double> set the thickness in cm of the concrete overburden\n";
    print "          (default is 300 and 0 means no overburden)\n";
    print "       -w <1> geometry with wires, <0> geometry with no wires\n";
    print "       -v <1> geometry with CRT, <0> geometry with no CRT \n";
    print "       -h prints this message, then quits\n";
    print "Remind: the file without wires has to be <filename_nowires.gdml> \n";
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
   <position name="PMTtube" unit="cm" x="0" y="0" z="@{[1.5*$PMTtube_x]}"/>
   <position name="posActiveInTPC"   unit="cm" x="$posTPCActive_x" y="$posTPCActive_y" z="$posTPCActive_z"/>
   <position name="posTPC0inCryo"    unit="cm" x="$TPCinCryo_x[0]" y="$TPCinCryo_y"    z="$TPCinCryo_z" />
   <position name="posCathode"  unit="cm" x="$posCat_x"  y="$TPCinCryo_y"    z="$TPCinCryo_z" />
   <position name="posTPC1inCryo"    unit="cm" x="$TPCinCryo_x[1]" y="$TPCinCryo_y"    z="$TPCinCryo_z" />
   <position name="posCryo1InWarmVessel"  unit="cm" x="$Cryo1InWarmVessel_x" y="$CryoInWarmVessel_y" z="$CryoInWarmVessel_z" />
   <position name="posCryo2InWarmVessel"  unit="cm" x="$Cryo2InWarmVessel_x" y="$CryoInWarmVessel_y" z="$CryoInWarmVessel_z" />
   <position name="posDetEncInWorld" unit="cm" x="$OriginXSet"     y="$DetEncl_yOffset"     z="$OriginZSet"/>
   <position name="posCenter"           unit="cm" x="0" y="0" z="0"/>
   <position name="posThermInsInDetEncl" unit="cm" x="$ThermInsInDetEncl_x" y="$ThermInsInDetEncl_y" z="$ThermInsInDetEncl_z"/>
   <position name="posWarmVesselInDetEncl" unit="cm" x="$WarmVesselInDetEncl_x" y="$WarmVesselInDetEncl_y" z="$WarmVesselInDetEncl_z"/>
   <position name="posCRTShellInDetEncl" unit="cm" x="0" y="$CRTSHELL_WV_OFFSET_Y" z="$CRTSHELL_WV_OFFSET_Z"/>
   <position name="posBuildingInWorld" unit="cm" x="0" y="@{[$Ground_y + $Building_y/2]}" z="0"/>
   <position name="posExpHallInWorld" unit="cm" x="0" y="@{[$Ground_y - $ExpHall_y/2 ]}" z="0"/>
   <rotation name="rPlus90AboutZPlus90AboutY"  unit="deg" x="0" y="90" z="90"/>
   <rotation name="rPlus90AboutX"       unit="deg" x="90" y="0" z="0"/>
   <rotation name="rMinus90AboutX"      unit="deg" x="-90" y="0" z="0"/>
   <rotation name="rPlus90AboutY"	unit="deg" x="0" y="90"   z="0"/>
   <rotation name="rPlus90AboutZ"	unit="deg" x="0" y="0"   z="90"/>
   <rotation name="rMinus90AboutZ"	unit="deg" x="0" y="0"   z="-90"/>
   <rotation name="rMinus90AboutY"      unit="deg" x="0" y="270" z="0"/>
   <rotation name="rMinus90AboutYMinus90AboutX"       unit="deg" x="270" y="270" z="0"/>
   <rotation name="rPlus90VAngleAboutX"	unit="deg" x="@{[90-$VAngle]}" y="0"   z="0"/>  
   <rotation name="rPlus90UAngleAboutX"	unit="deg" x="@{[90+$UAngle]}" y="0"   z="0"/>
   <rotation name="rPlus180AboutY"	unit="deg" x="0" y="180"   z="0"/>
   <rotation name="rIdentity"		unit="deg" x="0" y="0"   z="0"/>
   <rotation name="rPlusUAngleAboutX" 	unit="deg" x="$UAngle" y="0" z="0"/>
   <rotation name="rMinusVAngleAboutX" 	unit="deg" x="300" y="0" z="0"/>
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
  <element name="phosphorus" formula="P" Z="15"> <atom value="30.973"/>  </element>
  <element name="copper" formula="Cu" Z="29"> <atom value="63.5463"/>  </element>
  <element name="manganese" formula="Mn" Z="25"> <atom value="54.938043"/>  </element>
  <element name="vanadium" formula="V" Z="23"> <atom value="50.94151"/>  </element>

  <material name="Vacuum" formula="Vacuum">
   <D value="1.e-25" unit="g/cm3"/>
   <fraction n="1.0" ref="videRef"/>
  </material>

  <material name="ALUMINUM_Al" formula="ALUMINUM_Al">
   <D value="2.6990" unit="g/cm3"/>
   <fraction n="1.0000" ref="aluminum"/>
  </material>

  <material name="ALUMINUM_CRYO" formula="ALUMINUM_CRYO">
   <D value="0.5991" unit="g/cm3"/>
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

  <material name="STEEL_STAINLESS_Fe7Cr2Ni_WV" formula="STEEL_STAINLESS_Fe7Cr2Ni_WV">
   <!-- #<D value="0.3844" unit="g/cm3"/> -->
   <D value="0.38897" unit="g/cm3"/>
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

  <material name="Polystyrene">
   <D unit="g/cm3" value="1.06"/>
   <fraction n="0.077418" ref="hydrogen"/>
   <fraction n="0.922582" ref="carbon"/>
  </material>

  <material name="Polyurethane" formula="C27H36N2O10" >
   <D value="0.9" unit="g/cm3"/>
   <composite n="27" ref="carbon"/>
   <composite n="36" ref="hydrogen"/>
   <composite n="2" ref="nitrogen"/>
   <composite n="10" ref="oxygen"/>
  </material>

  <material name="STEEL_A992">
   <D unit="g/cm3" value="7.85"/>
   <fraction n="0.0022"  ref="carbon"/>
   <fraction n="0.005"   ref="copper"/>
   <fraction n="0.01"    ref="manganese"/>
   <fraction n="0.0044"  ref="nickel"/>
   <fraction n="0.00034" ref="phosphorus"/>
   <fraction n="0.0039"  ref="silicon"/>
   <fraction n="0.00044" ref="sulphur"/>
   <fraction n="0.001"   ref="vanadium"/>
   <fraction n="0.97272" ref="iron" />
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

# Set number of wires to default to zero, when $wires_on = 0, for a low memory 
# version. But if $wires_on = 1, calculate the number of wires on each side of each
# plane to be used in the for loops

my $NumberHorizontalWires = 0;
my $NumberCornerUWires = 0;
my $NumberCommonUWires = 0;
my $NumberCornerVWires = 0;
my $NumberCommonVWires = 0;


if ($wires_on == 1)
{
   # Number of wires in one corner
   #$NumberCornerVWires = int( $TPCWirePlane_y/$VWire_ypitch );
   #$NumberCornerWWires = int( $TPCWirePlane_y/$WWire_ypitch );
    $NumberCornerUWires = 480;
    $NumberCornerVWires = 480;


   # Total number of wires touching one vertical (longer) side
   # Note that the total number of wires per plane is this + another set of corner wires
   # $NumberSideUWires = int( $TPCWirePlane_z/$UWire_zpitch );
   # $NumberSideVWires = int( $TPCWirePlane_z/$VWire_zpitch );
   # $NumberSideWWires = int( $TPCWirePlane_z/$WWire_zpitch );

   # Number of wires per side that aren't cut off by the corner
   # $NumberCommonUWires = $NumberSideUWires - $NumberCornerUWires;
   # $NumberCommonVWires = $NumberSideVWires - $NumberCornerVWires;
    #$NumberCommonWWires = $NumberSideWWires - $NumberCornerWWires;
    $NumberCommonUWires  = 4640;
    $NumberCommonVWires  = 4640;

   # number of wires on the vertical plane
   #$NumberHorizontalWires = int( ($TPCWirePlane_y-$TPCWireThickness)/$UWirePitch );
   #Number of wires inthe Y plane-->Horizontal plane Induction I
    $NumberHorizontalWires = 1056;

}

# These XML comments throughout make the GDML file easier to navigate
print TPC <<EOF;

<!--+++++++++++++++++++ Y Wire Solids ++++++++++++++++++++++-->

EOF

if ($wires_on==1) 
{

#CommonWire = wires with same length

   print TPC <<EOF;
    <tube name="TPCWireYCommon"
      rmax="$TPCWireRadius"
      z="$TPCWirePlane_z"
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


<!--+++++++++++++++++++ U Wire Solids ++++++++++++++++++++++-->


EOF

# The corner wires for the U plane
if ($wires_on==1) 
{
#CORNER 
    $length = $CommonWireLength;

   for ($i = 0; $i < $NumberCornerUWires; ++$i)
    {
	 $length -= $DeltaLUCorner;

	print TPC <<EOF;

    <tube name="TPCWireU$i"
      rmax="$TPCWireRadius"
      z="$length"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>

EOF

#print(" $i $length \n");

   } #ends CORNER

   print TPC <<EOF;
    <tube name="TPCWireUCommon"
      rmax="$TPCWireRadius"
      z="$CommonWireLength"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

} else { 

print TPC <<EOF;

                   <!-- no wires in this GDML -->

EOF

}


print TPC <<EOF;


<!--+++++++++++++++++++ V Wire Solids ++++++++++++++++++++++-->


EOF

# The corner wires for the V plane
if ($wires_on==1) 
{
#CORNER
   $length = $CommonWireLength;

    for ($i = 0; $i < $NumberCornerVWires; ++$i)
   {

	$length -= $DeltaLVCorner;

   print TPC <<EOF;
    <tube name="TPCWireV$i"
      rmax="$TPCWireRadius"
      z="$length"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

#print(" $i $length \n");

    } #ends CORNER

   print TPC <<EOF;
    <tube name="TPCWireVCommon"
      rmax="$TPCWireRadius"
      z="$CommonWireLength"
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

  # Common Y wire logical volume, referenced many times
  print TPC <<EOF;
    <volume name="volTPCWireYCommon">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCWireYCommon" />
    </volume>
EOF

  # Corner U wires logical volumes
  for ($i = 0; $i < $NumberCornerUWires; ++$i)
  {
  print TPC <<EOF;
    <volume name="volTPCWireU$i">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCWireU$i" />
    </volume>
EOF

  } 

  # Common U wire logical volume, referenced many times
  print TPC <<EOF;
    <volume name="volTPCWireUCommon">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCWireUCommon" />
    </volume>
EOF

  # Corner V wires logical volumes
  for ($i = 0; $i < $NumberCornerVWires; ++$i)
  {
  print TPC <<EOF;
    <volume name="volTPCWireV$i">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCWireV$i" />
    </volume>
EOF

  }

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

#            ++++++++++++++++++++++  Y Plane  +++++++++++++++++++++++

# Create Y plane logical volume
print TPC <<EOF;


<!--+++++++++++++++++++++ Y Plane ++++++++++++++++++++++++-->


    <volume name="volTPCPlaneY">
      <materialref ref="LAr"/>
      <solidref ref="TPCPlane"/>
EOF

if ($wires_on==0)
{
print TPC <<EOF;

           <!-- no wires -->

EOF

} else {

    $ypos = ($NumberHorizontalWires+1)*$YWirePitch/2;

    for ($i = 0; $i < $NumberHorizontalWires; ++$i)
    {
	$ypos -= $YWirePitch;


print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireYCommon"/>
        <position name="posTPCWireY$i" unit="cm" x="0" y="$ypos " z="0"/>
        <rotationref ref="rIdentity"/>
      </physvol>
EOF

	#print("0 $ypos \n");

	#$ypos -= $YWirePitch;

    }

} #ends else


#            ++++++++++++++++++++++  U Plane  +++++++++++++++++++++++

# End U plane and create U plane logical volume
print TPC <<EOF;
    </volume>


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

#CORNERS

   $ypos1 = $TPCActive_y/2 ;
    $zpos1 = -$TPCActive_z/2 +  $CommonWireLength * $CosUAngle;

    $ypos2 = - $TPCActive_y/2 ;
    $zpos2 = -$TPCActive_z/2 ;

   # $ypos1 = $TPC_y/2 ;
   # $zpos1 = $TPC_z/2 ;

    #$ypos2 = -$TPC_y/2 ;
    #$zpos2 = $TPC_z/2 -$CommonWireLength * $CosUAngle;


   for ($i = 0; $i < $NumberCornerUWires; ++$i)
    {
	 $ypos1 += $UWire_ypitch  ;
	 $zpos2 -= $UWire_zpitch  ;

	 $ypos = ($ypos1+$ypos2)/2;
         $zpos = ($zpos1+$zpos2)/2;

	#print("U Corner wires: $i $zpos $ypos , ");
	

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireU$i"/>
        <position name="posTPCWireU$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotationref ref="rPlusUAngleAboutX"/>
      </physvol>
EOF

	$ypos = - $ypos;
	$zpos = - $zpos;


print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireU$i"/>
        <position name="posTPCWireU@{[$i+$NumberCommonUWires+$NumberCornerUWires]}" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotationref ref="rPlusUAngleAboutX"/>
      </physvol>
EOF

	#print("  $zpos $ypos \n");

    } #ends CORNER

#Common Wires

    #$zpos = -($NumberCommonUWires -1)*$UWire_zpitch/2;
   $zpos = (-$TPCActive_z +  $CommonWireLength * $CosUAngle) / 2.;
   #print("common wires $zpos \n");

    for ($i = 0; $i < $NumberCommonUWires; ++$i)
    {

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireUCommon"/>
        <position name="posTPCWireU$i" unit="cm" x="0" y="0 " z="$zpos"/>
        <rotationref ref="rPlusUAngleAboutX"/>
      </physvol>
EOF

    #print("U wires $i $zpos 0 \n");
	$zpos += $UWire_zpitch;

    }


} #ends else

#            ++++++++++++++++++++++  V Plane  +++++++++++++++++++++++

# End V plane and create V plane logical volume
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

#CORNERS

   $ypos1 = $TPCActive_y/2;
    $zpos1 = -$TPCActive_z/2;

    $ypos2 = -$TPCActive_y/2;
    $zpos2 = -$TPCActive_z/2 + $CommonWireLength * $CosVAngle;

   for ($i = 0; $i < $NumberCornerVWires; ++$i)
    {
	 $ypos1 -= $VWire_ypitch  ;
	 $zpos2 -= $VWire_zpitch  ;

	 $ypos = ($ypos1+$ypos2)/2;
         $zpos = ($zpos1+$zpos2)/2;

	#print("V Corner wires: $i $zpos $ypos , ");


print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireV$i"/>
        *<position name="posTPCWireV$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotationref ref="rMinusVAngleAboutX"/>
      </physvol>
EOF

	$ypos = - $ypos;
	$zpos = - $zpos;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireV$i"/>
        <position name="posTPCWireV@{[$i+$NumberCommonVWires+$NumberCornerVWires]}" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotationref ref="rMinusVAngleAboutX"/>
      </physvol>
EOF

	#print(" $zpos $ypos \n");

    } #ends CORNERS

#Common Wires

   #$zpos = -($NumberCommonVWires -1)*$VWire_zpitch/2;
   $zpos = (-$TPCActive_z +  $CommonWireLength * $CosVAngle) / 2.;

    for ($i = 0; $i < $NumberCommonVWires; ++$i)
    {

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireVCommon"/>
        <position name="posTPCWireV$i" unit="cm" x="0" y="0 " z="$zpos"/>
        <rotationref ref="rMinusVAngleAboutX"/>
      </physvol>
EOF

        #print("V wires: $i $zpos  0 \n");
	$zpos += $VWire_zpitch;

    }


} #ends else

print TPC <<EOF;
    </volume>
EOF

#+++++++++++++++++++++ Position physical wires Above +++++++++++++++++++++

my $VolY_x = (-$TPC_x/2) + 3*$WirePlaneSpacing; #+ $TPCWirePlane_x/2;    
my $VolU_x = (-$TPC_x/2) + 2*$WirePlaneSpacing; #+ $TPCWirePlane_x/2;    
my $VolV_x = (-$TPC_x/2) + 1*$WirePlaneSpacing; #+ $TPCWirePlane_x/2;  

#wrap up the TPC file
print TPC <<EOF;
    <volume name="volTPC">
      <materialref ref="LAr" />
      <solidref ref="TPC" />
     <physvol>
       <volumeref ref="volTPCPlaneY" />
       <position name="posTPCPlaneY" unit="cm" x="$VolY_x" y="0" z="0" />
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneU" />
       <position name="posTPCPlaneU" unit="cm" x="$VolU_x" y="0" z="0" />
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneV" />
       <position name="posTPCPlaneV" unit="cm" x="$VolV_x" y="0" z="0" />
     </physvol>
     <physvol>
       <volumeref ref="volTPCActive"/>
       <positionref ref="posActiveInTPC"/>
     </physvol>
     <physvol> 
       <volumeref ref="volRaceTrackTVolume"/>
       <positionref ref="posRaceTrackTInTPC"/>
       <rotationref ref="rIdentity"/>
     </physvol>
     <physvol> 
       <volumeref ref="volRaceTrackBVolume"/>
       <positionref ref="posRaceTrackBInTPC"/>
       <rotationref ref="rIdentity"/>
     </physvol>
     <physvol> 
       <volumeref ref="volRaceTrackUVolume"/>
       <positionref ref="posRaceTrackUInTPC"/>
       <rotationref ref="rPlus90AboutX"/>
     </physvol>
     <physvol> 
       <volumeref ref="volRaceTrackUVolume"/>
       <positionref ref="posRaceTrackDInTPC"/>
       <rotationref ref="rPlus90AboutX"/>
     </physvol>

    </volume>
</structure>
</gdml>
EOF

    close(GDML);

} #end of sub gen_TPC

##################################################################
#dg ##########  Race Tracks Parameters ###########################
##################################################################
sub define_RaceTrack()
{
# from "Design, construction and tests of the ICARUS T600 detector",
# ICARUS Collaboration (S. Amerio et al.). Jul 2004. 82 pp. 
# Nucl.Instrum.Meth. A527 (2004) 329-410
#
#  4 parts: Top (T) Bottom (B) Upstream (U) Downstream (D)
  $RT_epsilon = .0001 ;# extra safety space;
  $RaceTrack_d = 3.4;  # cm race track tube diameter
  $RaceTrack_TubeThick = 0.08; #cm race track tube thick
  $RaceTrack_lz = 1810.0; #cm length of race track structure
#  $RaceTrack_ly = 320.0 - 2*$RT_epsilon;  #cm length of race track tubes
  $RaceTrack_ly = 323.6201 ;  #cm height of race track structure (modified to be external to TPCActive_y
  $RaceTrack_number = 29; #number of race tracks
  $RaceTrack_pitch = 4.96; #distance between each race track tube
  $RaceTrack_ExternalRadius = $RaceTrack_d * 0.5;  # cm  external race track tube radius
  $RaceTrack_InnerRadius = $RaceTrack_d * 0.5-$RaceTrack_TubeThick;  # cm internal race track tube radius
  $RaceTrack_width = ($RaceTrack_pitch +$RT_epsilon) * $RaceTrack_number;
  $RaceTrackT_length = $RaceTrack_lz + $RT_epsilon; # guess about Horizontal length
  $RaceTrackTTube_length = $RaceTrackT_length - $RT_epsilon; # guess about Horizontal length
  $RaceTrackB_length = $RaceTrack_lz + $RT_epsilon; # guess about Horizontal length
  $RaceTrackBTube_length = $RaceTrackB_length - $RT_epsilon; # guess about Horizontal length
  $RaceTrackU_length  = $RaceTrack_ly - 2* $RaceTrack_d; # guess of vertical length
  $RaceTrackUTube_length  = $RaceTrackU_length - $RT_epsilon; # guess of vertical length
  $RaceTrack_thickness=$RaceTrack_d+2*$RT_epsilon; # 
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#dg +++++++++++++++++++++++++++++++++++ gen_RaceTracks +++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_RaceTracks{

#  4 parts: Top (T) Bottom (B) Upstream (U) Downstream (D)
# define the racetrack geometrical parameters
    define_RaceTrack();

    $RACETRACK= "icarus_racetracks" . $suffix . ".gdml";
    push (@gdmlFiles, $RACETRACK); # Add file to list of GDML fragments
    $RACETRACK = ">" . $RACETRACK;
    open(RACETRACK) or die("Could not open file $RACETRACK for writing");

# The standard XML prefix and starting the gdml
    print RACETRACK <<EOF;
<?xml version='1.0'?>
<gdml>
EOF

#Create the RACETRACK volumes:
# positions 

$RTy_T= ($RaceTrack_ly - $RaceTrack_thickness )/2. + $RT_epsilon ; 
$RTy_B=-($RaceTrack_ly - $RaceTrack_thickness )/2. + $RT_epsilon ;
$RTx_T= ($TPC_x - $RaceTrack_pitch -$RaceTrack_width )/2. ;
$RTx_B= $RTx_T ;
$RTx_U= $RTx_T ;
$RTz_U= -( $RaceTrackT_length - $RaceTrack_thickness ) / 2. - $RT_epsilon ;  
$RTx_D= $RTx_T ;
$RTz_D= -$RTz_U ;  
print RACETRACK <<EOF;
<define>

   <position name="posRaceTrackTInTPC" unit="cm" x="$RTx_T" y="$RTy_T" z="0"/>
   <position name="posRaceTrackBInTPC" unit="cm" x="$RTx_B" y="$RTy_B" z="0"/>
   <position name="posRaceTrackUInTPC" unit="cm" x="$RTx_U" y="0"  z="$RTz_U"/>
   <position name="posRaceTrackDInTPC" unit="cm" x="$RTx_D" y="0"  z="$RTz_D"/>

</define>   
<!--+++++++++++++++++++ RACETRACK Solids +++++++++++++++++++-->

EOF

print RACETRACK <<EOF;


<solids>
 <box name="RaceTrackTVolume"
  x="$RaceTrack_width"
  y="$RaceTrack_thickness"
  z="$RaceTrackT_length" 
  lunit="cm" />
 <box name="RaceTrackBVolume"
  x="$RaceTrack_width"
  y="$RaceTrack_thickness"
  z="$RaceTrackB_length" 
  lunit="cm"/>
 <box name="RaceTrackUVolume"
  x="$RaceTrack_width"
  y="$RaceTrack_thickness"
  z="$RaceTrackU_length" 
  lunit="cm"/>
 <tube name="RaceTrackTTubeVolume"
  rmax="$RaceTrack_ExternalRadius"
  rmin="$RaceTrack_InnerRadius"
  z="$RaceTrackTTube_length" 
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="RaceTrackBTubeVolume"
  rmax="$RaceTrack_ExternalRadius"
  rmin="$RaceTrack_InnerRadius"
  z="$RaceTrackBTube_length" 
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="RaceTrackUTubeVolume"
  rmax="$RaceTrack_ExternalRadius"
  rmin="$RaceTrack_InnerRadius"
  z="$RaceTrackUTube_length" 
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>

</solids>
EOF


# assume Upstream identical to Downstream part
# First define RaceTrack tubes 

print RACETRACK <<EOF;
<structure>

<!--+++++++++++++++++ RACETRACK Logical Volumes ++++++++++++++++++++-->

    <volume name="volRaceTrackTTubeVolume"> 
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>  
      <solidref ref="RaceTrackTTubeVolume"/>
    </volume>
    <volume name="volRaceTrackBTubeVolume"> 
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>  
      <solidref ref="RaceTrackBTubeVolume"/>
    </volume>
    <volume name="volRaceTrackUTubeVolume"> 
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>  
      <solidref ref="RaceTrackUTubeVolume"/>
    </volume>
EOF

# Define Upper part of RaceTrack
print RACETRACK <<EOF;

<!--+++++++++++++++++ RACETRACK Top Volume   ++++++++++++++++++++-->

    <volume name="volRaceTrackTVolume"> 
      <materialref ref="LAr"/>  
      <solidref ref="RaceTrackTVolume"/>


EOF


# Here it would be better to divide the volume into slices and then embedding tubes
# at the moment simply replicate each tube
$xTpos=-($RaceTrack_width - $RaceTrack_pitch - $RT_epsilon )/2.;
for ($it = 0; $it < $RaceTrack_number; $it++) 
{
  print RACETRACK <<EOF;
      <physvol>
        <volumeref ref="volRaceTrackTTubeVolume"/>
        <position name="posRTkTTube$it" unit="cm" x="$xTpos" y="0" z="0" />
        <rotationref ref="rIdentity" />
      </physvol>
EOF
  $xTpos+= $RaceTrack_pitch ;
}
print RACETRACK <<EOF;

    </volume>
EOF

# Define Lower part of RaceTrack
print RACETRACK <<EOF;

<!--+++++++++++++++++ RACETRACK Bottom Volume ++++++++++++++++++++-->

    <volume name="volRaceTrackBVolume"> 
      <materialref ref="LAr"/>  
      <solidref ref="RaceTrackBVolume"/>
EOF


# Here it would be better to divide the volume into slices and then embedding tubes
# at the moment simply replicate each tube
$xTpos=-($RaceTrack_width - $RaceTrack_pitch - $RT_epsilon )/2.;
for ($it = 0; $it < $RaceTrack_number; $it++) 
{
  print RACETRACK <<EOF;
      <physvol>
        <volumeref ref="volRaceTrackBTubeVolume"/>
        <position name="posRTkBTube$it" unit="cm" x="$xTpos" y="0" z="0" />
        <rotationref ref="rIdentity" />
      </physvol>
EOF
  $xTpos+= $RaceTrack_pitch ;
}
print RACETRACK <<EOF;

    </volume>
EOF

# Define Upstream part of RaceTrack
print RACETRACK <<EOF;

<!--+++++++++++++++++ RACETRACK Upstream Volume ++++++++++++++++++++-->

    <volume name="volRaceTrackUVolume"> 
      <materialref ref="LAr"/>  
      <solidref ref="RaceTrackUVolume"/>
EOF

# Here it would be better to divide the volume into slices and then embedding tubes
# at the moment simply replicate each tube
$xTpos=-($RaceTrack_width - $RaceTrack_pitch - $RT_epsilon )/2.;
for ($it = 0; $it < $RaceTrack_number; $it++) 
{
  print RACETRACK <<EOF;
      <physvol>
        <volumeref ref="volRaceTrackUTubeVolume"/>
        <position name="posRTkUTube$it" unit="cm" x="$xTpos" y="0" z="0" />
        <rotationref ref="rIdentity" />
      </physvol>
EOF
  $xTpos+= $RaceTrack_pitch ;
}
print RACETRACK <<EOF;

    </volume>
EOF

# close the RaceTrack structure 
print RACETRACK <<EOF;

<!--+++++++++++++++++ RACETRACK end structure +++++++++++++++++++-->
</structure>
EOF

#Close standard XML file
print RACETRACK <<EOF;
</gdml>
EOF

} 
#ends gen_RaceTracks

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_PMTs +++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_PMT {

    $PMT = "icarus_pmt" . $suffix . ".gdml";
    push (@gdmlFiles, $PMT); # Add file to list of GDML fragments
    $PMT = ">" . $PMT;
    open(PMT) or die("Could not open file $PMT for writing");

# The standard XML prefix and starting the gdml
    print PMT <<EOF;
<?xml version='1.0'?>
<gdml>
EOF

#Create the PMT volume original z=2.54
print PMT <<EOF;

<!--+++++++++++++++++++ PMT Solids ++++++++++++++++++++++-->

EOF
#<position name="PMTtube" unit="cm" x="30" y="0" z="0"/>
#was after print PMT <<EOF; line
#moved here to comment out for moving to position defs

print PMT <<EOF;
<solids> 
 <sphere name="PMTVolume"
  rmin="$PMTradiusInner"
  rmax="$PMTradiusOuter"
  deltaphi="360"
  deltatheta="90"
  aunit="deg"
  lunit="cm"/>
 <sphere name="PMTPassSphere"
  rmin="$PMTradiusInner"
  rmax="$PMTradiusOuter"
  deltaphi="360"
  deltatheta="90"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMTPassTube"
  rmax="$PMTtubeRmax"
  rmin="$PMTtubeRmin"
  z="$PMTtube_x" 
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <union name="PMTPassVolume">
  <first ref="PMTPassSphere"/>
  <second ref="PMTPassTube"/>
  <positionref ref="PMTtube"/>
  <rotationref ref="rMinus90AboutZ"/>
 </union>
</solids>   
EOF

#For some reasons, the Optical Sensitive Volume for PMT has to be LAr ... I found this info both in SBND and MicroBoone geometries
print PMT <<EOF;
<structure>

<!--+++++++++++++++++ PMT Logical Volumes ++++++++++++++++++++-->

    <volume name="volOpDetSensitive"> 
      <materialref ref="LAr"/>  
      <solidref ref="PMTVolume"/>
    </volume>
    <volume name="volNotOpDetSensitive"> 
      <materialref ref="Glass"/>  
      <solidref ref="PMTPassVolume"/>
    </volume>

</structure>
EOF

#Close standard XML file
print PMT <<EOF;
</gdml>
EOF

} 
#ends gen PMTs

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#am+++++++++++++++++++++++ gen_structure +++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Mech_Structure{

    $MECH = "icarus_Mech_Structure" . $suffix . ".gdml";
    push (@gdmlFiles, $MECH); # Add file to list of GDML fragments
    $MECH = ">" . $MECH;
    open(MECH) or die("Could not open file $MECH for writing");

# The standard XML prefix and starting the gdml
    print MECH <<EOF;
<?xml version='1.0'?>
<gdml>
EOF

#All the steel structure solids
print MECH <<EOF;

<!--+++++++++++++++++++ Structure Solids ++++++++++++++++++++++-->

EOF
print MECH <<EOF;
<solids> 
  <box name="LatExtMechBox"
  x="$ext_struct_x"
  y="$ext_struct_y"
  z="$ext_struct_z" 
  lunit="cm" />
  <box name="LatIntMechBox"
  x="@{[0.1+$int_struct_x]}"
  y="$int_struct_y"
  z="$int_struct_z" 
  lunit="cm" />
  <box name="Pillar"
  x="$int_struct_x"
  y="$int_struct_y"
  z="$int_struct_x" 
  lunit="cm" />
  <subtraction name="LatIntMechShell">
    <first ref="LatIntMechBox"/>
    <second ref="Pillar"/>
  </subtraction>
  <subtraction name="LatMechShell">
   <first ref="LatExtMechBox"/>
   <second ref="LatIntMechShell"/>
 </subtraction>
</solids>   
EOF

print MECH <<EOF;
<structure>

<!--+++++++++++++++++ Structure Logical Volumes ++++++++++++++++++++-->

<volume name="volLatMech"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>  
  <solidref ref="LatMechShell"/>
</volume>

</structure>
EOF

#Close standard XML file
print MECH <<EOF;
</gdml>
EOF

} 
#ends gen mechanical structure

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
      y="@{[$LAr_y + $GaseousAr_y]}"
      z="$TPC_z"/>
    <box name="GaseousArgon" lunit="cm" 
      x="$LAr_x"
      y="$GaseousAr_y"
      z="$LAr_z"/>
    <subtraction name="AlumShell">
      <first ref="Cryostat"/>
      <second ref="ArgonInterior"/>
    </subtraction>

    <box name="PMTPlane" lunit="cm" 
      x="$PMTPlane_x" 
      y="$PMTPlane_y" 
      z="$PMTPlane_z"/>

    <box name="Cathode" lunit="cm"
      x="$CPA_x"
      y="$CPA_y"
      z="$CPA_z"/>

    <box name="CathodeStrip" lunit="cm"
      x="$CPAStrip_x"
      y="$CPAStrip_y"
      z="$CPA_z"/>

</solids>
EOF

# Cryostat structure


print CRYO <<EOF;
<structure>
    <volume name="volAlumShell">
      <materialref ref="ALUMINUM_CRYO" />
      <solidref ref="AlumShell" />
    </volume>
    <volume name="volGaseousArgon">
      <materialref ref="ArGas"/>
      <solidref ref="GaseousArgon"/>
    </volume>
    <volume name="volCathodeStrip"> 
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>  
      <solidref ref="CathodeStrip"/>
    </volume>
    <volume name="volCathode">
      <materialref ref="LAr" />
      <solidref ref="Cathode" />
EOF

############################################################################################
#Positioning 64 cathode strips, 2.1 cm high.

$Num_CStrips = 64;
$yCSpos = -($CPA_y-$CPAStrip_pitch)/2;

    for ( $j=0; $j<$Num_CStrips; ++$j ){
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volCathodeStrip"/>
   <position name="posCStrip$j" unit="cm"  x="0" y="$yCSpos" z="0"/>
  </physvol>
EOF
$yCSpos+= $CPAStrip_pitch;
}
print CRYO <<EOF;
    </volume>
EOF
############################################################################################
print CRYO <<EOF;	
    <volume name="volPMTPlane">
      <materialref ref="LAr"/>
      <solidref ref="PMTPlane"/>
EOF

############################################################################################
#Positioning PMTs: positions from a file

#$PMT_x0 = 0; 
$PMT_x0 = $PMTradiusOuter/2 - $PMTGrid_offset/2 - $PMTtube_x ; 
$PMT_x1 = $PMT_x0 ;
@pmt_pos0 = read_pmt_pos("dispositionPMT.txt", $PMT_x0);
@pmt_pos1 = read_pmt_pos("dispositionPMT.txt", $PMT_x1);
$Num_PMTs0 = @pmt_pos0;

    for ( $i=0; $i<$Num_PMTs0; ++$i ){
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volOpDetSensitive"/>
   <position name="posPMT0$i" unit="cm" @pmt_pos0[$i]/>
   <rotationref ref="rPlus90AboutY"/>
  </physvol>
  <physvol>
   <volumeref ref="volNotOpDetSensitive"/>
   <position name="posPMT1$i" unit="cm" @pmt_pos1[$i]/>
   <rotationref ref="rMinus90AboutY"/>
  </physvol>
EOF
    }

############################################################################################
#Positioning Mechanical Structure elements:

#$TPCActive_z

$zMpos = -$TPCActive_z/2 + $ext_struct_z/2 ;
$mech_number = 9.; 
for ($im = 0; $im < $mech_number; $im++) 
{
  print CRYO <<EOF;
      <physvol>
        <volumeref ref="volLatMech"/>
        <position name="posLatMech$im" unit="cm" x="0" y="0" z="$zMpos" />
      <!--  <rotationref ref="rIdentity" /> -->
      </physvol>
EOF
  $zMpos+= $ext_struct_z ;
}
print CRYO <<EOF;
    </volume>
EOF

############################################################################################
print CRYO <<EOF;
    <volume name="volCryostat">
      <materialref ref="LAr" />
      <solidref ref="Cryostat" />
      <physvol>
        <volumeref ref="volGaseousArgon"/>
    <position name="posGaseousArgon" unit="cm" x="0" y="@{[$LAr_y/2]}" z="0" />
      </physvol>
      <physvol>
        <volumeref ref="volAlumShell"/>
      <position name="posAlumShell" unit="cm" x="0" y="0" z="0"/>
      </physvol>

     <physvol>
       <volumeref ref="volPMTPlane" />
       <position name="posPMTPlane1" unit="cm" x="@{[-$TPC_x - $PMTPlane_x/2 - $PMTWiresOffset ]}" y="$TPCinCryo_y" z="0" />
       <rotationref ref="rPlus180AboutY"/>
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

     <physvol>
       <volumeref ref="volPMTPlane" />
       <position name="posPMTPlane2" unit="cm" x="@{[$TPC_x + $PMTPlane_x/2 + $PMTWiresOffset ]}" y="$TPCinCryo_y" z="0" />
     </physvol>

EOF

print CRYO <<EOF;
    </volume>
</structure>
</gdml>
EOF

close(CRYO);
}


##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_CRT ++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_CRT()
{
# Create the CRT fragment file name,
# add file to list of output GDML fragments,
# and open it

# using a python script to generate the CRT geometry
    $CRTPYTHON="gen_crt_frags.py";
    my $ret=`python $CRTPYTHON `;
    

    $CRT = "icarus_crt" . $suffix . ".gdml";
    push (@gdmlFiles, $CRT);
# expecting the python script to generate a file "icarus_crt_test.gdml"
#    $CRT = ">" . $CRT;
#    open(CRT) or die("Could not open file $CRT for writing");


# The standard XML prefix and starting the gdml
#    print CRT <<EOF;
#<?xml version='1.0'?>
#<gdml>
#EOF


close(CRT);
}

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

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
     
    <box name="ThermIns" lunit="cm" 
      x="$ThermIns_x"
      y="$ThermIns_y"
      z="$ThermIns_z"/>

    <box name="WarmVessel" lunit="cm" 
      x="$WarmVessel_x"
      y="$WarmVessel_y"
      z="$WarmVessel_z"/>

    <box name="ThermInsInterior" lunit="cm" 
      x="@{[$ThermIns_x - 2*$FoamPadding]}"
      y="@{[$ThermIns_y - 2*$FoamPadding]}"
      z="@{[$ThermIns_z - 2*$FoamPadding]}"/>

    <box name="WarmVesselInterior" lunit="cm" 
      x="@{[$WarmVessel_x - 2*$WarmVesselThickness_x]}"
      y="@{[$WarmVessel_y - 2*$WarmVesselThickness]}"
      z="@{[$WarmVessel_z - 2*$WarmVesselThickness]}"/>

    <subtraction name="ThermInsShell">
      <first ref="ThermIns"/>
      <second ref="ThermInsInterior"/>
    </subtraction>

    <subtraction name="WarmVesselShell">
      <first ref="WarmVessel"/>
      <second ref="WarmVesselInterior"/>
    </subtraction>

</solids>
EOF


# Detector enclosure structure    
    print ENCL <<EOF;
<structure>

    <volume name="volThermIns">
      <materialref ref="Polyurethane"/> 
      <solidref ref="ThermInsShell"/>
    </volume>

    <volume name="volWarmVessel">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni_WV"/>
      <solidref ref="WarmVesselShell"/>
    </volume>

    <volume name="volDetEnclosure">
      <materialref ref="Air"/>
      <solidref ref="DetEnclosure"/>

    <physvol>
      <volumeref ref="volThermIns"/>
      <positionref ref="posThermInsInDetEncl"/>
    </physvol>

    <physvol>
      <volumeref ref="volCRT_Shell"/>
      <positionref ref="posCRTShellInDetEncl"/>
    </physvol>

    <physvol>
      <volumeref ref="volWarmVessel"/>
      <positionref ref="posWarmVesselInDetEncl"/>
    </physvol>

    <physvol>
      <volumeref ref="volCryostat"/>
      <positionref ref="posCryo1InWarmVessel"/>
    </physvol>

    <physvol>
      <volumeref ref="volCryostat"/>
      <positionref ref="posCryo2InWarmVessel"/>
    </physvol>
	
EOF
if ($crt_on==1) 
{
#+++++++++++++++++++++++++ Position CRT modules ++++++++++++++++++++++++++
    #$YPosMINOSSide = $Cryostat_y/2+$MINOSModWidth/2;


} #end of crt_on 

 else { 

print ENCL <<EOF;

        <!-- This GDML version has no crt and uses much less memory -->

EOF

}

print ENCL <<EOF;

</volume>

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

#Dirt global
#    <box name="DirtWorld" lunit="cm" 
#      x="$World_x" 
#      y="@{[$World_y/2 + $Ground_y]}" 
#      z="$World_z"/>

# All the World solids.
print WORLD <<EOF;
<solids>
    <box name="World" lunit="cm" 
      x="$World_x" 
      y="$World_y" 
      z="$World_z"/>

    <box name="BoxDirtAlongX" lunit="cm" 
      x="@{[$Hall_x-0.1]}" 
      y="@{[$ExpHall_y-0.1]}" 
      z="@{[$World_z/2-$Hall_z/2-0.1]}"/>

    <box name="BoxDirtAlongZ" lunit="cm" 
      x="@{[$World_x/2-$Hall_x/2-0.1]}" 
      y="@{[$ExpHall_y-0.1]}" 
      z="@{[$World_z-0.1]}"/>

    <box name="BoxDirtBottom" lunit="cm" 
      x="@{[$World_x-0.1]}" 
      y="@{[$World_y/2+$Ground_y-$ExpHall_y]}" 
      z="@{[$World_z-0.1]}"/>

EOF

if ($concrete_on != 0)
{ 

#part of the overburden above the dirt
print WORLD <<EOF;
    <box name="OverburdenUp" lunit="cm" 
      x="$Overburden_x" 
      y="@{[2*$Overburden_y/3]}" 
      z="$Overburden_z"/>

    <box name="OverburdenDown" lunit="cm" 
      x="$Overburden_x" 
      y="@{[$Overburden_y/3]}" 
      z="$Overburden_z"/>
EOF
}

else {
print WORLD <<EOF;
        <!-- This GDML version has no overburden-->
EOF
}

print WORLD <<EOF;

    <box name="Building" lunit="cm" 
      x="$Hall_x" 
      y="$Building_y" 
      z="$Hall_z"/>

    <box name="AirBuilding" lunit="cm" 
      x="@{[$Hall_x - $HallWallThicnekss]}" 
      y="@{[$Building_y - $HallWallThicnekss/2]} " 
      z="@{[$Hall_z - $HallWallThicnekss]}"/>

    <subtraction name="WallBuilding">
     <first ref="Building"/>
     <second ref="AirBuilding"/>
     <position name="posAirBuilding" unit="cm" x="0" y="@{[-$HallWallThicnekss/4]}" z="0"/>
     </subtraction>

    <box name="ExpHall" lunit="cm" 
      x="$Hall_x" 
      y="$ExpHall_y" 
      z="$Hall_z"/>

    <box name="AirExpHall" lunit="cm" 
      x="@{[$Hall_x - $HallWallThicnekss]}" 
      y="@{[$ExpHall_y - $HallWallThicnekss/2]}" 
      z="@{[$Hall_z - $HallWallThicnekss]}"/>

    <subtraction name="WallExpHall">
     <first ref="ExpHall"/>
     <second ref="AirExpHall"/>
     <position name="posAirExpHall" unit="cm" x="0" y="@{[$HallWallThicnekss/4]}" z="0"/>
     </subtraction>
</solids>
EOF


# World structure: Building + underground experimental hall + detector

print WORLD <<EOF;
<structure>
EOF

#Building: building + upper part of the overburden

if ( $concrete_on != 0) { 
print WORLD <<EOF;
    <volume name="volOverburdenUp" >
      <materialref ref="Air"/>
      <solidref ref="OverburdenUp"/>
    </volume>
EOF
}

print WORLD <<EOF;  

    <volume name="volBuilding" >
      <materialref ref="Air"/>
      <solidref ref="AirBuilding"/>
    </volume>

    <volume name="volWallBuilding" >
      <materialref ref="Concrete"/>
      <solidref ref="WallBuilding"/>

EOF

#if ($concrete_on !=0) {

#print WORLD <<EOF;    

#     <physvol>
#        <volumeref ref="volOverburdenUp"/>
#	<position name="posOverburdenUp" unit="cm" x="$posOverburden_x" y="@{[-$Building_y/2 + $Overburden_y/6 ]}" z="$posOverburden_z" />
#      </physvol>
#EOF

#}

#Experimental hall: experimental hall + overburden

print WORLD <<EOF; 
</volume>
    <volume name="volExpHall" >
      <materialref ref="Air"/>
      <solidref ref="AirExpHall"/>
    </volume>

EOF

if ( $concrete_on != 0) { 
print WORLD <<EOF;
    <volume name="volOverburdenDown" >
      <materialref ref="Air"/>
      <solidref ref="OverburdenDown"/>
    </volume>
EOF
}

print WORLD <<EOF;

    <volume name="volWallExpHall" >
      <materialref ref="Concrete"/>
      <solidref ref="WallExpHall"/>
EOF


#if ( $concrete_on != 0) { 
#print WORLD <<EOF;
#     <physvol>
#        <volumeref ref="volOverburdenDown"/>
#	<position name="posOverburdenDown" unit="cm" x="$posOverburden_x" y="@{[$ExpHall_y/2 - $Overburden_y/3]}" z="$posOverburden_z" />
#      </physvol>
#EOF
#}



#Complete world: building + underground parts + detector

print WORLD <<EOF;
</volume>

    <volume name="volBoxDirtAlongX" >
      <materialref ref="Dirt"/>
      <solidref ref="BoxDirtAlongX"/>
    </volume>

    <volume name="volBoxDirtAlongZ" >
      <materialref ref="Dirt"/>
      <solidref ref="BoxDirtAlongZ"/>
    </volume>

    <volume name="volBoxDirtBottom" >
      <materialref ref="Dirt"/>
      <solidref ref="BoxDirtBottom"/>
     </volume>

    <volume name="volWorld" >
      <materialref ref="Air"/>
      <solidref ref="World"/>

    <physvol>
    <volumeref ref="volWallBuilding"/>
    <positionref ref="posBuildingInWorld"/>
    </physvol>
    
    <physvol>
    <volumeref ref="volWallExpHall"/>
    <positionref ref="posExpHallInWorld"/>
    </physvol>
    
     <physvol>
        <volumeref ref="volOverburdenUp"/>
	<position name="posOverburdenUp" unit="cm" x="$posOverburden_x" y="@{[$Ground_y + $Overburden_y/3 ]}" z="$posOverburden_z" />
      </physvol>

     <physvol>
        <volumeref ref="volOverburdenDown"/>
	<position name="posOverburdenDown" unit="cm" x="$posOverburden_x" y="@{[$Ground_y - $Overburden_y/6]}" z="$posOverburden_z" />
      </physvol>

      <physvol>
        <volumeref ref="volBoxDirtAlongX"/>
	<position name="posBoxDirtAlongX1" unit="cm" x="$posOverburden_x" y="@{[$Ground_y-0.5*($ExpHall_y-0.1)]}" z="@{[$World_z/2-0.5*($World_z/2-$Hall_z/2)]}"/>
      </physvol>

      <physvol>
        <volumeref ref="volBoxDirtAlongX"/>
	<position name="posBoxDirtAlongX2" unit="cm" x="$posOverburden_x" y="@{[$Ground_y-0.5*($ExpHall_y-0.1)]}" z="@{[-$World_z/2+0.5*($World_z/2-$Hall_z/2)]}"/>
      </physvol>

      <physvol>
        <volumeref ref="volBoxDirtAlongZ"/>
	<position name="posBoxDirtAlongZ1" unit="cm" x="@{[$World_x/2- 0.5*($World_x/2-$Hall_x/2)]}" y="@{[$Ground_y-0.5*($ExpHall_y-0.1)]}" z="$posOverburden_z"/>
      </physvol>

      <physvol>
        <volumeref ref="volBoxDirtAlongZ"/>
	<position name="posBoxDirtAlongZ2" unit="cm" x="@{[-$World_x/2+0.5*($World_x/2-$Hall_x/2)]}" y="@{[$Ground_y-0.5*($ExpHall_y-0.1)]}" z="$posOverburden_z"/>
      </physvol>

      <physvol>        
        <volumeref ref="volBoxDirtBottom"/>
	<position name="posBoxDirtBottom" unit="cm" x="$posOverburden_x" y="@{[-$World_y/2+0.5*($World_y/2+$Ground_y-$ExpHall_y)]}" z="$posOverburden_z"/>
      </physvol>     

      <physvol>
        <volumeref ref="volDetEnclosure"/>
	<positionref ref="posDetEncInWorld"/>
      </physvol>

    </volume>
</structure>
</gdml>
EOF

#big dirt volume
#<physvol>
#<volumeref ref="volDirtWorld"/>
#<position name="posDirtWorld" unit="cm" x="$posOverburden_x" y="@{[-$World_y/4 + $Ground_y/2]}" z="$posOverburden_z"/>
#</physvol>

# make_gdml.pl will take care of <setup/>

close(WORLD);
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++ write_fragments ++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#sub gen_CRT()
#{
# Create the TPC fragment file name,
# add file to list of output GDML fragments,
# and open it
#    $CRT = "crt.gdml";
#    push (@gdmlFiles, $CRT);
#}

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

