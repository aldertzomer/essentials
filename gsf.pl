#!/usr/local/bin/perl

use strict ;
use warnings ;

use fgweb::general::general ;
#use fgweb::general::cmd_input ;
use fgweb::var ;

# global parameters: DO NOT change
# =====================================

#my %ini ;
#$ini{'-i'}{desc} = 'Input file.' ;
#$ini{'-i'}{checkisexisting} = 1 ;
#$ini{'-i'}{checkisstring} = 1 ;
#$ini{'-i'}{mandatory} = 1 ;

my $vars = var::fgweb_vars() ;
my $tool = $vars->{cgi}{root}.'gsf/gsf_run.pl' ;

# functions
# =====================================

sub parseparams {       # shows the usage information if applicable
#    my $title = "\nBy Sacha van Hijum [s.a.f.t.van.hijum\@run.nl]\n" ;
#    cmd_input::check_cmd_input(\%ini, \@ARGV,$title) ;
} # parseparams

sub runtool {
    my $shortname = $tool ;
    $shortname = $1 if ($shortname =~ /.+\/(.)/) ;

    $ENV{'PATHPROGS'} = $vars->{cgi}{root}.'gsf/' ;
    $ENV{'TOUCHLOC'} = $vars->{sys}{touch} ; 
    $ENV{'ZIP_LOC'} = $vars->{sys}{zip} ;
    $ENV{'R'} = $vars->{tools}{R} ;
    $ENV{'PS2PNG'} = $vars->{tools}{ps2png} ;
    
    system($tool) == 0 or general::error('errors during run of '.$shortname) ;
} # runtool

# Main
# =====================================

parseparams() ;
runtool() ;
