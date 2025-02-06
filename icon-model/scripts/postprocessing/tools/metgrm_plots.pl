#!/usr/bin/perl

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#
=head1 NAME

B<metgrm_plots.pl> - Plot ICON meteograms with C<ncl>.

=head1 SYNOPSIS

B<metgrm_plots.pl> [B<-c> colormap] [B<-D>] B<-e> var1[,var2,...] [B<-i>] [B<-I> ID]
        [B<-l> lev0] [B<--oTypye> outputType] [B<-s> station,[station2,...]]
        [B<-S>] [B<-t> tRange] [B<-T> title] [B<-z> zAxis] F<meteogram_file>

=head1 DESCRIPTION

Plot ICON meteograms with ncl-scripts B<mtgrm_plot.ncl> or B<mtgrm_plot_sfc.ncl>.

=head1 OPTIONS

=over 5

=item B<-c> colormap

NCL colormap for contour plots. Default: C<BlueDarkRed18>.
Use the name C<Red2Blue> to reverse it.
C<WhiteBlue> is color map C<BlueWhite> reversed.

=item B<-D>  or  B<--debug>

Debug option. Repeat for higher debug level. B<ncl> is not executed for debug levels
3 or higher.

=item B<-e> var[,var2,var3,...]

Plot variables var, var2, var3, ...
Use
C<var=all> to plot all variables,
C<var=3d> or C<var=var_name> to plot all 3D-variables,
C<var=2d> or C<var=sfcvar_name> to plot all 2D-variables,
C<var=tile> to plot all tile variables, or
C<var=notile> to plot all 2D-variables except for tile variables.
If C<var=/rex/> or C<var=/rex/i> then all variables matching the regular expression
C<rex> are plotted.

=item B<-h>

This help.

=item B<-i>

List variables and stations in the meteogram file.

=item B<-I> ID

ID for plots.

=item B<-l> lev0

Plot 3D-variables starting at level C<lev0>.
Negativ numbers means to plot levels 0 to nlev-1+lev0
Default: 0.

=item B<--oType> outputType

Output type.  Possible values: eps, png, ps
Default: C<png>.

=item B<-s> station[,station2,...]

Make plot for these stations.

=item B<-S>

Only list the stations in the meteogram file.

=item B<-t> tRange

Plot only the given time range. C<tRange=date0:date1> where
date0 and date1 have the format C<YYYYmmddHH> or only C<YYYYmmdd>
or the ICON date format C<YYYYmmddTHH>. Leaving out date0 or date1
defaults to start or end of the date range in the meteogram.
C<tRange=date> is equivalent to C<tRange=:date>.

=item B<-T> title

Title for the plot. The default title is C<ICON>.

=item B<-z> zAxis

Use C<zAxis> as vertical coordinate for contour plots. Default is C<ml> 
for model levels. Possible other values C<h> for height.

=back

=head2 Note

=head2 created by:  Helmut P. Frank, DWD, FE 13,  28.01.2019

=head2  Modifications:

=over

=item 16.04.2021

Use the ncl-scripts in directory $NCL_SCRIPT_DIR to make the plots. The default is
C<NCL_SCRIPT_DIR=.> .

=item 31.01.2019

Add option B<-t>.

=back

=cut

use Getopt::Long;
Getopt::Long::Configure( qw(bundling permute) );
use strict;

my $prog = $0;
$prog =~ s?^.*/??;

my $debug;
my $help;
my $id;
my $info;
my @vars;
my @stations;
my $list_stations;
my $oType = 'png';
my $tRange;
my $title = 'ICON';
my $colormap = 'BlueDarkRed18';
my $lev0 = 0;
my $zAxis = 'ml';

my $result = GetOptions( "colormap|c=s" => \$colormap,
                         "Debug|D+"   => \$debug,
                         "e|variable=s" => \@vars,
                         "help|h"      => \$help,
                         "info|i"      => \$info,
                         "lev0|l=i"    => \$lev0,
                         "oType=s"     => \$oType,
                         "ID|I=s"      => \$id,
                         "station|s=s" => \@stations,
                         "S"           => \$list_stations,
                         "tRange|t=s"  => \$tRange,
                         "title|T=s"   => \$title,
                         "zAxis|z=s"   => \$zAxis)
    or die "\tError calling GetOptions in $prog\n";

my $metgram_file = $ARGV[0];

@stations = split( /,/, join( ',', @stations));
@vars     = split( /,/, join( ',', @vars));
$tRange = check_tRange( $tRange);

print STDERR "Options:\n",
             "\tdebug     = $debug\n",
             "\tvariables = @vars\n",
             "\thelp      = $help\n",
             "\tinfo       = $info\n",
             "\tID        = $id\n",
             "\tstation   = @stations\n",
             "\tlist_stations= $list_stations\n",
             "\ttRange      = $tRange\n",
             "\ttitle       = $title\n",
             "\tmetgrm_file = $metgram_file\n",
    if ($debug);

unless ( defined($metgram_file)) {
    print STDERR "No meteogram file specified!\n";
    $help = 1;
}

exec ("perldoc $0") if ($help);

die "Meteogram file $metgram_file not readable!\n" unless ( -r $metgram_file);

open NC, "ncdump -v station_name,var_name,sfcvar_name,date $metgram_file |" or
    die "Error opening pipe from ncdump -v station_name,var_name,sfcvar_name,date $metgram_file";

my @station_name = read_nc_variable( 'station_name', 'char', $debug);
my @var_name     = read_nc_variable( 'var_name', 'char', $debug);
my @sfcvar_name  = read_nc_variable( 'sfcvar_name', 'char', $debug);
my @date         = read_nc_variable( 'date', 'char', $debug);
close NC;

if ( $info) {
    print "ICON meteogram file: $metgram_file\n";
    write_info( "Stations", 3, sort @station_name);
    write_info( "Variables", 4, sort @var_name);
    write_info( "Surface variables", 4, sort @sfcvar_name);
    exit;
}
if ( $list_stations) {
    print "ICON meteogram file: $metgram_file\n";
    write_info( "Stations", 3, sort @station_name);
    exit;
}

@stations = @station_name if ( $stations[0] eq 'all');

#  Find the variables to plot
my $var;
my @var_plot;
foreach $var (@vars) {
    VAR: {
        if ( $var eq 'all') {
            push @var_plot, ( @var_name, @sfcvar_name);
	    next VAR;
        }
        if ( $var eq '3d' || $var eq 'var_name') {
            push @var_plot, @var_name;
	    next VAR;
        }
        if ( $var eq '2d' || $var eq 'sfcvar_name') {
            push @var_plot, @sfcvar_name;
	    next VAR;
        }
        if ( $var =~ /^tile/) {
            push ( @var_plot, ( grep( /_T_/, @sfcvar_name)) );
	    next VAR;
        }
        if ( $var =~ /^notile/) {
            push ( @var_plot, ( grep( !/_T_/, @sfcvar_name)) );
	    next VAR;
        }
        if ( $var =~ /^\w+$/) {
            push @var_plot, $var;
        } elsif ( $var =~ s[^/(.+)/$][$1]) {
            push @var_plot, ( grep( /$var/, @vars, @sfcvar_name));
        } elsif ( $var =~ s[^/(.+)/i$][$1]) {
            push @var_plot, ( grep( /$var/i, @vars, @sfcvar_name));
        }
    }
}
my %var_plot;
foreach $var (@var_plot) {
    $var_plot{$var} += 1;
}
@var_plot = keys %var_plot;

my %station_no;
my %var_name;
my %sfcvar_name;
my ( $x, $i);
$i = 1;
foreach $x (@station_name) {
    $station_no{$x} = $i++;
}
$i = 1;
foreach $x (@var_name) {
    $var_name{$x} = $i++;
}
$i = 1;
foreach $x (@sfcvar_name) {
    $sfcvar_name{$x} = $i++;
}

unless ( defined($id)) {
    my $cymdg = substr( $date[0],0, 11);
    $cymdg =~ s/T//;
    $id = "icon$cymdg";
}

$title = $metgram_file unless ( defined($title));

my ( $oFile, $ncl_arg, $ncl_script, $ncl_cmd, $station);
$title =~ s/ /\\ /g;

my $ncl_dir = defined( $ENV{NCL_SCRIPT_DIR}) ? $ENV{NCL_SCRIPT_DIR} : '.';
foreach $station ( @stations) {
    unless ( defined($station_no{$station})) {
	print STDERR "Station $station is not in meteogram file $metgram_file!\n";
	next;
    } elsif ( $debug) {
	print STDERR "iStation=$station_no{$station} for $station\n";
    }
    foreach $var (@var_plot) {
	$oFile = $id . "_$station" . "_$var.$oType";
	$ncl_arg = qq(iFile="$metgram_file" oFile="$oFile" oType="$oType" varName="$var" iStation=$station_no{$station} expnum="$title") .
                   qq( levMin="$lev0" zAxis="$zAxis" tRange="$tRange" colormap="$colormap");
	$ncl_arg =~ s/"/\\"/g;

	if ( defined($var_name{$var})) {
            $ncl_script = $ncl_dir . "/mtgrm_plot.ncl"
	} elsif ( defined($sfcvar_name{$var})) {
            $ncl_script = $ncl_dir . "/mtgrm_plot_sfc.ncl"
	} else {
            print STDERR "Variable $var not in meteogram file $metgram_file!\n";
	    next;
        }
        die "ncl-script $ncl_script not readable!" unless ( -r $ncl_script);
        print STDERR qq(ncl -n $ncl_script $ncl_arg\n) if ($debug);
	if ( $debug <= 2) {
            $ncl_cmd = qq(ncl -n $ncl_script $ncl_arg);
            $result = 0xffff & system( $ncl_cmd);
	    if ( $result == 0xff00) {
		print STDERR "ncl failed! Error executing\n",
		print STDERR "$ncl_cmd\n";
            } elsif ( $result > 0x80) {
                $result >>= 8;
                print STDERR "ncl ran with non-zero exit status $result\n" unless ($result==1);
            } elsif ($result & 0x80) {
                $result &= ~0x80;
                print STDERR "ncl ran with coredump from signal $result\n"
	    } else{
 		print "Created file $oFile\n\n";
	    }
	}
    }
}

sub read_nc_variable {
#
#   Read netCDF variable $variable of $type from netcdf file which was opened
#   at file handle NC.
#
    my ( $variable, $type, $debug, $missing) = @_;
    my ( $i, $x);
    my @values;
    while (<NC>) {
        next if (/^\s*$/);
        if ( s/$variable\s*=//) {
            chomp( $x = $_);
            while ( $x !~ /;\s*$/) {
                print STDERR if ($debug>2);
                chomp($_ = <NC>);
                $x .= $_;
            }
            if ( $type eq 'char') {
                @values = split( m{"\s*,\s*"}, $x);
                $values[0]  =~ s{^\s*"}{};
                $values[-1] =~ s{"\s*;\s*$}{};
            } else {
                $x =~ s/[,;]/ /g;
		$x =~ s/$missing/?/g if ( defined($missing));
                @values = split( ' ', $x);
            }
            print STDERR "$variable (", scalar(@values), "): @values\n" if ($debug>1);
            last;
        }
    }
    return @values;
}

sub write_info {
    my ( $text, $m, @list) = @_;
    print "$text:\n";
    my $i = 0;
    my $l;
    while ( defined( $l = shift @list)) {
	next unless ( length($l));
        print "\t $l";
	print "\n" if ( ++$i%$m == 0);
    }
    print "\n";
}

sub check_tRange {
#
#  Convert from ICON data and time format or from only dates to format YYYYmmddHH
#
    my ( $t_range, $dummy) = @_;
    my @t_range = split( /[\s:]/, $t_range);
    foreach my $tr ( @t_range) {
        if ( $tr =~ /\*/ || $tr !~ /\S/ ) {
            $tr = '';
        } elsif ( $tr =~ /(\d{8})T?(\d+)?/i) {
            my ( $d, $t) = ( $1, $2);
	    if ( length($t) >= 2) {
                $tr = $d . substr($t,0,2);
            } else {
                $tr = $d . sprintf( "%2.2d", $t);
            }
        } else {
            die  "Invalid date and time $tr!\n";
        }
    }
    if ( $#t_range >= 1) {
	$t_range =  "$t_range[0]:$t_range[1]";
    } else {
        if ( $t_range !~ /:/) {
	    $t_range =  $t_range[0];
        } elsif ( $t_range =~ /^\s:/) {
	    $t_range =  ":$t_range[0]";
        } else {
	    $t_range =  "$t_range[0]:";
        }
    }
    $t_range = '' if ( $t_range eq ':');
    return $t_range;
}
