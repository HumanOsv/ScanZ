#!/usr/bin/perl -s

use strict;
use warnings;
use Benchmark;

# Large in ángstrom 
my $largeAxisZ = 60; 
my $window     = 0.5;

# Mopac Config
my $semiempirical = "AUX LARGE PM6-D3H4 NSPA=78.4";
my $Charge = 0;
my $Multiplicity = 1;
my $ncpus = 4;
my $mem = "8GB";


# # #
# Path for software mopac 
my $path_bin_mopac = "/opt/mopac/MOPAC2016.exe";
#
#
my %Atomic_number = ( '89'  => 'Ac', '13'  => 'Al', '95'  => 'Am', '51'  => 'Sb',	
	                  '18'  => 'Ar', '33'  => 'As', '85'  => 'At', '16'  => 'S',  
					  '56'  => 'Ba', '4'   => 'Be', '97'  => 'Bk', '83'  => 'Bi',	
                      '107' => 'Bh', '5'   => 'B', 	'35'  => 'Br', '48'  => 'Cd',	
	                  '20'  => 'Ca', '98'  => 'Cf',	'6'   => 'C',  '58'  => 'Ce',	
	                  '55'  => 'Cs', '17'  => 'Cl',	'27'  => 'Co', '29'  => 'Cu',	
	                  '24'  => 'Cr', '96'  => 'Cm', '110' => 'Ds', '66'  => 'Dy',
	                  '105' => 'Db', '99'  => 'Es', '68'  => 'Er', '21'  => 'Sc',	
	                  '50'  => 'Sn', '38'  => 'Sr', '63'  => 'Eu', '100' => 'Fm',	
	                  '9'   => 'F',  '15'  => 'P',  '87'  => 'Fr', '64'  => 'Gd',	
	                  '31'  => 'Ga', '32'  => 'Ge', '72'  => 'Hf', '108' => 'Hs',	
                      '2'   => 'He', '1'   => 'H',  '26'  => 'Fe', '67'  => 'Ho',	
					  '49'  => 'In', '53'  => 'I',  '77'  => 'Ir', '70'  => 'Yb',
					  '39'  => 'Y',  '36'  => 'Kr', '57'  => 'La', '103' => 'Lr',	
					  '3'   => 'Li', '71'  => 'Lu', '12'  => 'Mg', '25'  => 'Mn',	
                      '109' => 'Mt', '101' => 'Md', '80'  => 'Hg', '42'  => 'Mo',	
					  '60'  => 'Nd', '10'  => 'Ne', '93'  => 'Np', '41'  => 'Nb',	
					  '28'  => 'Ni', '7'   => 'N',  '102' => 'No', '79'  => 'Au',	
					  '76'  => 'Os', '8'   => 'O', 	'46'  => 'Pd', '47'  => 'Ag',	
					  '78'  => 'Pt', '82'  => 'Pb',	'94'  => 'Pu', '84'  => 'Po',	
					  '19'  => 'K',  '59'  => 'Pr', '61'  => 'Pm', '91'  => 'Pa',	
					  '88'  => 'Ra', '86'  => 'Rn', '75'  => 'Re', '45'  => 'Rh',	
					  '37'  => 'Rb', '44'  => 'Ru', '104' => 'Rf', '62'  => 'Sm',
					  '106' => 'Sg', '34'  => 'Se', '14'  => 'Si', '11'  => 'Na',
					  '81'  => 'Tl', '73'  => 'Ta', '43'  => 'Tc', '52'  => 'Te',	
					  '65'  => 'Tb', '22'  => 'Ti', '90'  => 'Th', '69'  => 'Tm',	
					  '112' => 'Uub','116' => 'Uuh','111' => 'Uuu','118' => 'Uuo',	
					  '115' => 'Uup','114' => 'Uuq','117' => 'Uus','113' => 'Uut',
					  '92'  => 'U',  '23'  => 'V',  '74'  => 'W',  '54'  => 'Xe',
                      '30'  => 'Zn', '40'  => 'Zr' );



sub energy_per_mol {
	# filename
	my ($input_file) = @_;
	#
	my @HeaderLines;
	open(HEADER,"$input_file.arc") or die "Unable to open $input_file.arc";
	@HeaderLines  = <HEADER>;
	close HEADER;
	#
	my $energy_m;
	while (my $HLine = shift (@HeaderLines)) {
		chomp ($HLine);
		my $hatfield_1 = "TOTAL ENERGY";
		if ( $HLine =~/$hatfield_1/ ){
			my $energy  = $HLine;
			my @words_1 = split (" ",$energy);
			$energy_m   = $words_1[3];
		}
	}
	return $energy_m;
}
###################################
# Read files
sub read_file    {
	# filename
	my ($input_file) = @_;
	my @array        = ();
	# open file
	open(FILE, "<", $input_file ) || die "Can't open $input_file: $!";
	while (my $row = <FILE>) {
		chomp($row);
		push (@array,$row);
	}
	close (FILE);
	# return array	
	return @array;
}
###################################
# Reference file
sub format_xyz   {
	my ($input_file) = @_;
	#
	my @array_coord = ();
	#
	my $tam = scalar (@{$input_file});
	for ( my $i = 2 ; $i < $tam ; $i = $i + 1 ){
		if ( length(@$input_file[$i]) > 2) { 
			my @array_tabs  = split (/\s+/,@$input_file[$i]);
			my $radii_val;
			if ( exists $Atomic_number{$array_tabs[0]} ) {
				# exists
				$radii_val = $Atomic_number{$array_tabs[0]};
			} else {
				# not exists
				$radii_val = $array_tabs[0] ;
			}
			my $strong = "$radii_val\t$array_tabs[1]\t$array_tabs[2]\t$array_tabs[3]";
			push (@array_coord,$strong);
		}
	}
	return @array_coord;
}
###################################
# Reaction coordinate 
sub z_axis_path  {
    my $cell_center = 0; 
    my $total_number_of_cell_z = ( $largeAxisZ / $window );
    $total_number_of_cell_z=~ s/\.\d+$//;
    #
    my @array = ();
    #
    my $div_points_z = ($total_number_of_cell_z / 2);
    $div_points_z=~ s/\.\d+$//;
    my $neg_points_z = $div_points_z * -1;
    my $pos_points_z = $div_points_z;
    #
    for (my $z=$neg_points_z; $z <= ($pos_points_z + 0.0001); $z++) { 
        my $coord_z = ( ( $z * $window ) + $cell_center );
        my $c_z = sprintf '%.6f', $coord_z;
        push (@array,$c_z);
    }
    return @array;
}
###################################
# Input file Mopac
sub MopacInput   {
	#
	my $filebase     = $_[0];
	my $coordsMM     = $_[1];
	my $MopacInput   = "$filebase.mop";
	my $iteration    = $_[2];
	my $Headerfile   = $_[3];
	my $Charge       = $_[4];
	my $Multiplicity = $_[5];
	#
	my $mem          = $_[7];
	#
	my $tmp   = "+1";
	my @words = split (/\n/,$coordsMM);
	#
	open (COMFILE, ">$MopacInput");
	#
	my $word;
	# Spin multiplicity:
	if ( $Multiplicity == 0 ) { $word = "NONET"   };			
	# singlet	- 0 unpaired electrons
	if ( $Multiplicity == 1 ) { $word = "SINGLET" };
	# doublet	- 1 unpaired electrons
	if ( $Multiplicity == 2 ) { $word = "DOUBLET" };
	# triplet	- 2 unpaired electrons
	if ( $Multiplicity == 3 ) { $word = "TRIPLET" };
	# quartet	- 3 unpaired electrons
	if ( $Multiplicity == 4 ) { $word = "QUARTET" };
	# quintet	- 4 unpaired electrons			
	if ( $Multiplicity == 5 ) { $word = "QUINTET" };
	# sextet	- 5 unpaired electrons
	if ( $Multiplicity == 6 ) { $word = "SEXTET"  };
	# septet	- 6 unpaired electrons
	if ( $Multiplicity == 7 ) { $word = "SEPTET"  };
	# octet	- 7 unpaired electrons
	if ( $Multiplicity == 8 ) { $word = "OCTET"   };
	#
	my $ncpus        = ($_[6] * 2);
	if ( $ncpus == 0 ) {
		print COMFILE "$Headerfile $word CHARGE=$Charge";
	} else {
		# The maximum number of threads is normally equal to the number of cores, 
		# even if each core supports two threads.
		# In the special case of THREADS=1, parallelization is switched off.
		print COMFILE "$Headerfile $word CHARGE=$Charge THREADS=$ncpus";
	}	
	print COMFILE "\n";
	print COMFILE "Scan job $iteration\n";
	print COMFILE "\n";	
	foreach my $i (@words){
		my @axis    = split (" ",$i);
		#
		my $label  = $axis[0];  
		my $axis_x = $axis[1];
		my $axis_y = $axis[2];
		my $axis_z = $axis[3];
		#
		print COMFILE "$label\t$axis_x\t$tmp\t$axis_y\t$tmp\t$axis_z\t$tmp\n";
	}
	print COMFILE "\n";
	print COMFILE "\n";
	close (COMFILE);
	#
	return $MopacInput;
}
###################################
# Energy ouputs mopac
sub energy_mopac {
	my ($num_atoms_xyz, $name_file) = @_;	
	# directorio
	my $dir = './';
	#
	my @array = ();
	my @array_coords_mopac = ();
	# abrir directorio
	opendir(DIR, $dir) or die $!;
	#
	while (my $file = readdir(DIR)) {
		# Use a regular expression to ignore files beginning with a period
		next if ($file =~ m/^\./);
		next unless ($file =~ m/\.arc$/);
		# sin extension
		my $fileSnExt = $file;
		$fileSnExt =~ s/\..*$//;
		#
		if ( -e "$fileSnExt.arc" ) {
			push( @array, "$fileSnExt.arc");
		} else {
			print "WARNING: Mopac file $file error termination\n";		
		}
		#
	}
	closedir(DIR);
	#
	my $tam_esc = scalar (@array);
	if ($tam_esc == 0) { print "ERROR problem MOPAC $tam_esc files .arc, check .out\n"; exit(0);}
	#
	my @HeaderLines = ();
	my @ZeroPoint   = ();
	my @energyy     = ();
	my $energy      = '';
	my $number_atoms = $num_atoms_xyz;
	foreach my $i (@array) {
		#
		open(HEADER,"$i") or die "Unable to open $i";
		@HeaderLines  = <HEADER>;
		close HEADER;
		#
		while (my $HLine = shift (@HeaderLines)) {
			chomp ($HLine);
			#
			my $hatfield_1 = "TOTAL ENERGY";
			if ( $HLine =~/$hatfield_1/ ){
				$energy = $HLine;
				my @words_1 = split (" ",$energy);
				push (@energyy,$words_1[3]); 
			}
		}
	}
	#
	my @idx = sort { $array[$a] cmp $array[$b] } 0 .. $#array;
	my @energyy_1     = @energyy[@idx];
	my @array_1       = @array[@idx];
	#
	my $mopac_all_geome = $name_file; 
	open (ARCFILE, ">$mopac_all_geome");
	#
	for ( my $i = 0; $i < scalar(@array_1); $i++) {
		# Convert kcal/mol
		#my $kcal         = 23.0605419453;
		#my $theBest_E    = abs($energyy_1[0] - $energyy_1[$i]);
		# Total Energy
		#my $convert_E    = ($theBest_E * $kcal);	
		#
		open (HEADER,"$array_1[$i]") or die "Unable to open $i";
		my @HeaderLines = <HEADER>;
		close HEADER;
		#
		my $count_lines = 0;
		my $first_line  = 0;
		my @array_lines = ();
		while (my $HLine = shift (@HeaderLines)) {
			my $hatfield = "FINAL GEOMETRY OBTAINED";
			if ( $HLine =~/$hatfield/ ){
				$first_line = $count_lines;			
			}
			$count_lines++;
			push (@array_lines,$HLine);
		}
		my $tmp_rest = $first_line + 3;	
		my $lala     =  $count_lines;
		#
		my $concat;
		#
		print ARCFILE "$number_atoms\n";
		print ARCFILE "$energyy_1[$i] eV ";
		#my $number = sprintf '%05f', $convert_E;
		#print ARCFILE "$number Kcal/mol\n";
		#
		print ARCFILE "$array_1[$i] \n";
		#
		for ( my $i = $tmp_rest; $i < $count_lines; $i++) {
			my $wordlength = length ($array_lines[$i]);    
			#
			if ( $wordlength > 3) {
				chomp ($array_lines[$i]);
				my @words = split (" ",$array_lines[$i]);
				print ARCFILE "$words[0]\t$words[1]\t$words[3]\t$words[5]\n";
				$concat.= "$words[0]\t$words[1]\t$words[3]\t$words[5]\n";
			}
		}
		push (@array_coords_mopac,$concat);
	}
	close ARCFILE;
	return @energyy_1;
}
#
#
my ($file_name1,$file_name2) = @ARGV;
if (not defined $file_name1) {
	die "\nScanMopac must be run with:\n\nUsage:\n\tScanMopac.pl [XYZ-Mol1] [XYZ-Mol2]\n\n\n";
	exit(1);
}
#
if (not defined $file_name2) {
	die "\nScanMopac must be run with:\n\nUsage:\n\tScanMopac.pl [XYZ-Mol1] [XYZ-Mol2]\n\n\n";
	exit(1);
}
# Initial Time
my $tiempo_inicial  = new Benchmark;
#
# Read and parse files1
my @data_file1    = read_file($file_name1);
my @coords_file1  = format_xyz (\@data_file1);
#
my $without_extension;
($without_extension = $file_name1) =~ s/\.[^.]+$//;
# Read and parse files2
my @data_file2    = read_file($file_name2);
my @coords_file2  = format_xyz (\@data_file2);
#
my $without_extension_2;
($without_extension_2 = $file_name2) =~ s/\.[^.]+$//;
# Calculo de energia por cada molecula
my $str_1;
my $str_2;
foreach my $i (@coords_file1) {
	my ($Atom,$axisX,$axisY,$axisZ)  = split (" ",$i);
	$str_1.="$Atom\t$axisX\t$axisY\t$axisZ\n";
}
foreach my $j (@coords_file2) {
	my ($Atom,$axisX,$axisY,$axisZ)  = split (" ",$j);
	$str_2.="$Atom\t$axisX\t$axisY\t$axisZ\n";
}
#
my $MopacInput_1 = MopacInput ($without_extension,$str_1,$without_extension,$semiempirical,$Charge,$Multiplicity,$ncpus,$mem);
my $MopacInput_2 = MopacInput ($without_extension_2,$str_2,$without_extension_2,$semiempirical,$Charge,$Multiplicity,$ncpus,$mem);
#
system ("$path_bin_mopac $MopacInput_1 >tmp_mopac_1.txt 2>tmp_mopac_2.txt");
system ("$path_bin_mopac $MopacInput_2 >tmp_mopac_1.txt 2>tmp_mopac_2.txt");
#
my $Emol1 = energy_per_mol ($without_extension);
my $Emol2 = energy_per_mol ($without_extension_2);
#
# Delete tmp files Part 1
system ("rm *.out *.mop *.aux");
system ("rm *.arc");
system ("rm tmp_mopac_1.txt 2>tmp_mopac_2.txt");
#
my @array_Z     = z_axis_path();
# 
# Total Number of Atoms
my $system1Numb = $data_file1[0];
my $system2Numb = $data_file2[0];
my $NumberAtoms = $system1Numb + $system2Numb;
# Animation
#
my @FinalCoords = ();
for (my $z=0; $z < scalar (@array_Z); $z++) {
    my $stringTmp;
    #print "$NumberAtoms\n";
    #print "\n";
    foreach my $i (@coords_file1) {
        my ($Atom,$axisX,$axisY,$axisZ)  = split (" ",$i);
        my $sumZ = $axisZ + $array_Z[$z];
        #print "$Atom\t$axisX\t$axisY\t$sumZ\n";
        $stringTmp.="$Atom\t$axisX\t$axisY\t$sumZ\n";
    }
    foreach my $i (@coords_file2) {
        my ($Atom,$axisX,$axisY,$axisZ)  = split (" ",$i);
        #print "$Atom\t$axisX\t$axisY\t$axisZ\n";
        $stringTmp.="$Atom\t$axisX\t$axisY\t$axisZ\n";
    }
    push(@FinalCoords,$stringTmp);
}
#
my $uniq_count   = 0;
for (my $i=0; $i < scalar (@array_Z); $i++) {
    my $id            = sprintf '%.6d', $uniq_count;
    my $filebaseMopac = "Scan$id";
    my $MopacInput = MopacInput ($filebaseMopac,$FinalCoords[$i],$i,$semiempirical,$Charge,$Multiplicity,$ncpus,$mem);
    #
    system ("$path_bin_mopac $MopacInput >tmp_mopac_1.txt 2>tmp_mopac_2.txt");
    #
    $uniq_count++;
}
print "Energy\n";
my  @Energy_mopac = energy_mopac ($NumberAtoms,"Coords-$without_extension.xyz");
#
#my @energy_sort = sort  { $a <=> $b } @Energy_mopac;
#
for (my $i=0; $i < scalar (@array_Z); $i++) {
#	my $kcalmol = (abs ($energy_sort[0]) - abs ($Energy_mopac[$i]) ) * 23.0605419453 ;
	my $kcalmol = ($Energy_mopac[$i] - ($Emol1 + $Emol2)) * 23.0605419453;
	my $number = sprintf '%05f', $kcalmol;
	#
	print "$number\t$array_Z[$i]\n";
}
my $tiempo_final = new Benchmark;
my $tiempo_total = timediff($tiempo_final, $tiempo_inicial);
print "\n\tExecution Time: ",timestr($tiempo_total),"\n";
print "\n";

# Delete tmp files Part 2
system ("rm *.out *.mop *.aux");
system ("rm *.arc");
system ("rm tmp_mopac_1.txt 2>tmp_mopac_2.txt");

