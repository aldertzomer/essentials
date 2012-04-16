package add2db_link ;

use fgweb::var ;

my $vars ;
my @options = ();
$vars = var::fgweb_vars() ;

my $add2db_name = 'add2db' ;
my $add2db_url  = $vars->{url}{root}.'add2db/add2db_start.php' ;
my $add2db_logo = $vars->{url}{root}.'add2db/images/add2db.gif' ;
my $add2db_display_name = 'Add2DB' ;

sub add2db_html_link
{       
    my $referer_tool    = shift() ;
    my $output_html     = shift() ;
    my $output_html_url = shift() ;
    my $file            = shift() ;
    my $filetype        = shift() ;
    my $format          = shift() ;
    my $datatype        = shift() ;
    @options = ();
    
    push (@options,  shift() );
    push (@options,  shift() );
    push (@options,  shift() );
    push (@options,  shift() );
    push (@options,  shift() );
    
                        
    my $curoption = undef;
    my $success = 1 ;
    
    if ( (defined $referer_tool) and (defined $output_html_url) and (defined $output_html) and (defined $file) and (defined $filetype) ) {
	# ok, we have the necessary variables
    } else {
	$success = 0 ;
    }
    if ($success == 1) {
	open (OUT,'>'.$output_html) or $success = 0 ;
    }
    if ($success == 1) {
	print OUT '<html>'."\n";
	print OUT '<body>'."\n";
	print OUT '<form name="autosend" enctype="multipart/form-data" action="'.$add2db_url.'" method="post" />'."\n";
        print OUT '<input type="hidden" name="'.$add2db_name.'_'.$filetype.'_select" value="yes" />'."\n";
        print OUT '<input type="hidden" name="referer_tool" value="'.$referer_tool.'" />'."\n";
        print OUT '<input type="hidden" name="'.$add2db_name.'_'.$filetype.'_fname" value="'.$file.'" />'."\n";
	if (defined $format) {
    	    print OUT '<input type="hidden" name="'.$add2db_name.'_'.$filetype.'_format" value="'.$format.'" />'."\n";
	}
	if (defined $datatype) {
    	    print OUT '<input type="hidden" name="'.$add2db_name.'_'.$filetype.'_type" value="'.$datatype.'" />'."\n";
	}
	if (defined @options) {
	    foreach $curoption (@options) {
		if ($curoption =~ /^(.+)\=(.+)/) {
    		    print OUT '<input type="hidden" name="'.$add2db_name.'_'.$filetype.'_'.$1.'" value="'.$2.'" />'."\n";
		}
	    }
	}
        print OUT '</form>'."\n";
        print OUT '<script language="javascript">'."\n";
        print OUT '<!--'."\n";
        print OUT 'document.forms["autosend"].submit();'."\n";
        print OUT '-->'."\n";
        print OUT '</script>'."\n";
        print OUT '</body>'."\n";
        print OUT '</html>'."\n";
        close(OUT) ;
    }
    if ($success == 1) {
	$out  = '<a href="'.$output_html_url.'" target="myadd2db" onclick=\'window.open( "", "myadd2db", "width=820,height=500,resizable,scrollbars,screenX=20,screenY=40,left=20,top=40");\' >'."\n" ;
	$out .= '<img src="'.$add2db_logo.'" border="0" alt="'.$add2db_display_name.'" /></a>'."\n" ;
    } else {
	$out = '' ;
    }
    return $out ;
}


1;
