#!/usr/local/bin/perl
# GSF perl script. Calls scripts for parsing the genbank (module 1)  and generating the insertion sites (module 2), 
# parses configfile, downloads reads, does alignments  (module 3)  and runs statistical analysis module 4 (R)
# prerequisites: R, GNU tools, curl, pass. fq_all2std_noqual.pl, fastx_barcode_splitter_trimmer, gbk2pttfna.pl, tasite.pl
# fastx_barcode_splitter_trimmer requires Text::LevenshteinXS if you want to split your barcodes quickly

use strict;
use warnings;
use fgweb::general::general;
use fgweb::general::add2db_link;
use fgweb::var;

# global variables

#inputfile from website
my $configfile = 'configfile.txt';
my $genbank = $ENV{'TEMPLATE'};
my $pwd = $ENV{'PWD'};
my $insertion = $ENV{'INSERTIONSITE'};
my $sessiondir = $ENV{'SESSION_DIR'};
my $workdir = $ENV{'WORK_DIR'};
my $normalization = $ENV{'NORMALIZATION'};
my $barcodemismatch = $ENV{'BARCODEMISMATCH'};
my $minhit = $ENV{'MINHIT'};
my $strand = $ENV{'STRAND'};
my $librarysize = 2 * $ENV{'LIBRARYSIZE'};
my $statmethod = $ENV{'STATMETHOD'};


#fgweb
my $vars  = var::fgweb_vars() ;
my $Rnew  = $vars -> {tools}{Rnew} ;
my $pass  = $vars -> {tools}{pass} ;
my $curl  = $vars -> {sys}{curl} ;
my $cat   = $vars -> {sys}{cat} ;
my $cut   = $vars -> {sys}{cut} ;
my $mv    = $vars -> {sys}{mv} ;
my $grep  = $vars -> {sys}{'grep'} ;
my $sed   = $vars -> {sys}{sed} ;
my $tr    = $vars -> {sys}{'tr'} ;
my $sort  = $vars -> {sys}{'sort'} ;
my $uniq  = $vars -> {sys}{uniq} ;
my $paste = $vars -> {sys}{paste} ;
my $perl  = $vars -> {tools}{perl} ;
#my $perl = '/usr/local/bin/perl';

#perl scripts for parsing the genbank, the insertion site finder, the fastx tool kit and maq

my $gbk2fna 		= $vars->{cgi}{root}.'gsf/gbk2pttfna.pl';
my $tasitescript 	= $vars->{cgi}{root}.'gsf/tasite.pl';
my $fq_all2std  	= $vars->{cgi}{root}.'gsf/fq_all2std_noqual.pl';
my $fastx_split 	= $vars->{cgi}{root}.'gsf/fastx_barcode_splitter_trimmer.pl';

# R scripts for stats

my $gene_edger 		= $vars->{cgi}{root}."gsf/R/".$statmethod."_gene_edgeR.R";
my $esgene_edger 	= $vars->{cgi}{root}."gsf/R/".$statmethod."_essentialgenes_edgeR.R";
my $ta_edger 		= $vars->{cgi}{root}."gsf/R/".$statmethod."_ta_edgeR.R";


######################################################################################
# Run programs down here. Comment out sections to speed up a rerun or to debug.

# Parsing genbank
print "Parsing selected genbank file. \n \n";
qx("$perl" "$gbk2fna" -gbk="$genbank" -s=genome.fasta -sa=genome.ptt 2>> gbk2fna.error);

# Finding insertion sites
print "Searching TA sites in fasta file of selected genome.\n \n";
qx("$perl" "$tasitescript" genome.fasta "$insertion" >genome.ta 2>> tasitescript.error);

# Run the analysis
get_and_convert() ;

# Run the statistics
doStats();

# Generate html
generatehtml() ;

# debug option
qx(chmod 777 -R "$workdir");

######################################################################################
# get_and_convert
# Main subroutine. Everything else is started from here. 
# This subroutine will parse the configfile, download the read files in fastq or
# fasta format and will convert all files to fasta. Also it will generate the required 
# files for other subroutines

sub get_and_convert {
    my @downloadedfiles;
    my @samples;
    my $barcode;
    my $inputformat;
    my $inputlink;
    my $filecount;
    my $inputlinkprevious;
    my $lineconfig;
    my $transposon;
    my $sampletype;
    my $samplecount;
    my $samples;
    my $inputcompression;
    my $decompress;
    my $library ;
    my %config ;
    my @splot ;
    my @splut ;
    #unused    my $transposonmismatch;

    print "Preparing headers of output files\n";
    writeheaders();
    
    # parse the configfile line by line. each line contains the link, the barcode, the transposon sequence, the sampletype (control or target)  and the input format
    # If the barcode or transposon is not present, use N as sequence. 
    
    print "Parsing configfile\n";
    open(CONFIGFILE, $configfile) or general::error("Can't open file:" .$configfile);
    $filecount = 0;
    $samplecount = 0;
    $inputlinkprevious = "a";
    while (<CONFIGFILE>){
        $lineconfig = $_ ;
        chomp $lineconfig ;
        @splot = split '\t', $lineconfig;
        $inputlink = $splot[0] ;
	$barcode = $splot[1] ;
	$transposon = $splot[2] ;
	$sampletype = $splot[3] ;
	$library =  $splot[4] ;
	$inputformat = $splot[5];
	$inputcompression = $splot[6];
        next if ($inputlink eq "Link");
        next if ($inputlink eq "link");
        next if ($inputlink eq "#");
        next if ($inputlink eq "#link");
        next if ($inputlink eq "#Link");
        $samplecount = $samplecount + 1;
        #check if we have not downloaded this file before
        if ($inputlink ne $inputlinkprevious) {
	    $filecount = $filecount +1;
	
	    #define readfile and sample
	    my $readfile = "reads".$filecount.".fasta";
	    my $sample = "split_".$sampletype."_sample".$samplecount."";

	    #Sanity checks
	    if ($barcode =~ /[^ATCGNatcgn]/) { general::error("barcode ".$barcode." contains non ATCGN characters!\n");}
	    if ($transposon =~ /[^ATCGNatcgn]/) {general::error("Transposon sequence ".$transposon." contains non ATCGN characters!\n");}
	    
	    #Decompression options. SRA has bz2 compressed files and now we can directly pull stuff of SRA.
	    if ($inputcompression eq "none") {$decompress = " "};
	    if ($inputcompression eq "bz2")  {$decompress = "| bzcat "};
	    if ($inputcompression eq "gz") {$decompress = "| zcat "};
	    if ($inputcompression eq "zip")  {$decompress = "| funzip "};
	    
	    #Download and convert read files	     
	    print "Downloading $inputlink to $readfile, converting it from $inputformat to fasta\n";
	    if ($inputformat eq "fasta") { qx("$curl" "$inputlink" --stderr "$readfile".error $decompress > "$readfile");}
	    if ($inputformat eq "fastq") { qx("$curl" "$inputlink" --stderr "$readfile".error $decompress |"$fq_all2std" fq2fa   	> "$readfile");}
	    if ($inputformat eq "scarf") { qx("$curl" "$inputlink" --stderr "$readfile".error $decompress |"$fq_all2std" scarf2std	|"$fq_all2std" fq2fa  > "$readfile");}
	    if ($inputformat eq "export"){ qx("$curl" "$inputlink" --stderr "$readfile".error $decompress |"$fq_all2std" export2std	|"$fq_all2std" fq2fa  > "$readfile");}
	    if ($inputformat eq "solid") { qx("$curl" "$inputlink" --stderr "$readfile".error $decompress |"$fq_all2std" csfa2std	|"$fq_all2std" fq2fa  > "$readfile");}

	    general::error("download of ".$inputlink." failed\n") if (-z $readfile);
	    general::error("conversion of ".$inputlink." to fasta failed\n") unless (-e $readfile);
	    
	    #create barcode file for later use
    	    if ($barcode ne "N") {
    		print "File $readfile has barcode(s) assigned. First sample is $sample\n";
    		my $barcodefile = "barcodes".$filecount.".txt";
    		open (BARCODE, "> $barcodefile");
        	print BARCODE "".$sampletype."_sample".$samplecount."\t".$barcode."\n";
        	close (BARCODE);
        	push (@downloadedfiles, $filecount);
            }
            
            #if there's no barcode, there's also no need to split. Immediately give the readfile a sample name instead of a number                                  
	    if ($barcode eq "N") {
		    print "File $readfile has no barcode assigned. Assuming the entire file is sample $sample\n";
		    qx("$mv" "$readfile" transposon_"$sample".fasta 2>> get_and_convert_mv.error);
	    }	    
	    
	    #create transposon site file for later use (long or short, has an effect on downstream processing)
	    if (length $transposon < 8) {
    		if ($transposon ne "N") {
    		    my $transposonfile = "short_transposon_".$sample.".txt";
		    open (TRANSPOSON, "> $transposonfile");
		    print TRANSPOSON "$sample\t$transposon\n";
		    close (TRANSPOSON);
		}
	    }
	    if (length $transposon >= 8) {
    		my $transposonfile = "long_transposon_".$sample.".fasta";
		open (TRANSPOSON, "> $transposonfile");
		print TRANSPOSON "\>transposon\n";
		print TRANSPOSON "$transposon\n";
		close (TRANSPOSON);
	    }
	    
	    #add sample to file with control and target samples for gene edgeR
	    my $countsgene = "gene_".$sample.".counts.txt";
	    open (TARGETSGENES, ">> gene_targets.txt");
	    print TARGETSGENES "".$countsgene."\t".$sampletype."\t".$sampletype."\t".$library."\n";
	    close (TARGETSGENES);
	    
	    #add sample to file with all samples as control samples for essential genes edgeR
	    open (TARGETSESSENTIALGENES, ">> essentialgenes_targets.txt");
	    print TARGETSESSENTIALGENES "".$countsgene."\ttarget\ttarget\t".$library."\n";
	    close (TARGETSESSENTIALGENES);
	    
	    #add sample to file with control and target samples for ta site edgeR
	    my $tacounts = "ta_".$sample.".counts.txt";
	    open (TARGETSTA, ">> ta_targets.txt");
	    print TARGETSTA "".$tacounts."\t".$sampletype."\t".$sampletype."\t".$library."\n";
	    close (TARGETSTA);
            push(@samples, $sample);
	    $inputlinkprevious = $inputlink;
	}

	# If the readfile from the inputlink contains another sample, process it similar as above, but do some extra sanity checks
	
	else {
	    
	    #define readfiles and samples
	    my $readfile = "reads".$filecount.".fasta";
	    my $sample = "split_".$sampletype."_sample".$samplecount."";
	    print "$readfile also contains sample $sample\n";
	                
	    #Sanity checks
	    if ($barcode =~ /[^ATCGNatcgn]/) {general::error("Barcode ".$barcode." contains non ATCGN characters\n");}
	    if ($barcode eq "N")  {general::error("Your previous sample for file ".$readfile." did contain a barcode while this one does not. Check your config files! \n");}
	    unless (-e $readfile) {general::error("Your previous sample for file ".$readfile." did not contain a barcode while this one does. Check your config files! \n");}
	    if ($transposon =~ /[^ATCGNatcgn]/) {general::error("Transposon sequence ".$transposon." contains non ATCGN characters\n");}
	    
	    #add barcode to barcode file
	    my $barcodefile = "barcodes".$filecount.".txt";
	    open (BARCODE, ">> $barcodefile");
	    print BARCODE "".$sampletype."_sample".$samplecount."\t".$barcode."\n";
	    close (BARCODE);
	    
	    #create transposon site file for later use (long or short, has an effect on downstream processing)
	    if (length $transposon < 8) {
    		if ($transposon ne "N")  {
    		    my $transposonfile = "short_transposon_".$sample.".txt";
		    open (TRANSPOSON, "> $transposonfile");
		    print TRANSPOSON "$sample\t$transposon\n";
		    close (TRANSPOSON);
		}
	    }
	    if (length $transposon >= 8) {
    		my $transposonfile = "long_transposon_".$sample.".fasta";
		open (TRANSPOSON, "> $transposonfile");
		print TRANSPOSON "\>transposon\n";
		print TRANSPOSON "$transposon\n";
		close (TRANSPOSON);
	    }

	    #add sample to file with control and target samples for gene edgeR
	    my $countsgene = "gene_".$sample.".counts.txt";
	    open (TARGETSGENES, ">> gene_targets.txt");
	    print TARGETSGENES "".$countsgene."\t".$sampletype."\t".$sampletype."\t".$library."\n";
	    close (TARGETSGENES);
	    
	    #add sample to file with all samples as control samples for essential genes edgeR
	    open (TARGETSESSENTIALGENES, ">> essentialgenes_targets.txt");
	    print TARGETSESSENTIALGENES "".$countsgene."\ttarget\ttarget\t".$library."\n";
	    close (TARGETSESSENTIALGENES);
	    
	    #add sample to file with control and target samples for ta site edgeR
	    my $tacounts = "ta_".$sample.".counts.txt";
	    open (TARGETSTA, ">> ta_targets.txt");
	    print TARGETSTA "".$tacounts."\t".$sampletype."\t".$sampletype."\t".$library."\n";
	    close (TARGETSTA);
            push(@samples, $sample);
	}
    }
    #enter the numbers of the downloaded files into the trimming and splitting routine creating per sample a fasta file
    trim_and_split(@downloadedfiles);

    #remove the long or short transposon if necessary
    remove_transposon(@samples);

    # align the fasta files of the samples to the ta sites, count the number of reads per ta site and count the number of reads per gene
    pass_alignments(@samples);

    #calculate number of unique TA sites per gene
    doTAcounts();
    
    #annotate the TA sites
    doTAannotation();

}#END get_and_convert

######################################################################################
# writeheaders
# this subroutine prepares some initial files.

sub writeheaders {
    #create initial gene_targets.txt readDGE object input file
    open (TARGETSGENE, '>gene_targets.txt');
    print TARGETSGENE "files \t group \t description \t library \n";
    close (TARGETSGENE);
    #create initial essentialgenes_targets.txt readDGE object input file
    open (TARGETSESSENTIALGENES, '>essentialgenes_targets.txt');
    print TARGETSESSENTIALGENES "files \t group \t description \t library \n";
    close (TARGETSESSENTIALGENES);
    #create initial ta_targets.txt readDGE object input file
    open (TARGETSTA, '>ta_targets.txt');
    print TARGETSTA "files \t group \t description \t library \n";
    close (TARGETSTA);
}# END writeheaders

#################################################################################
# trim_and_split
# This subroutine will split the file on barcodes. Replace with other splitter if 
# needed. Output must be named transposon_$sample.fasta

sub trim_and_split(@) {
    
    my @dlarray=@_;
    my $dlarrayref;

    foreach $dlarrayref (@dlarray) {
    	print "Splitting reads".$dlarrayref.".fasta\n";
    	stat("barcodes".$dlarrayref.".txt");
    	if (-e _) {
    	    qx("$cat" reads"$dlarrayref".fasta | $fastx_split --bcfile barcodes"$dlarrayref".txt --prefix transposon_split_ --suffix .fasta --bol --mismatches "$barcodemismatch"  2>>trim_and_split_fastxsplit.error);
	    unlink("reads".$dlarrayref.".fasta");
	}
	else {
	    #just rename the read file to the sample name. Normally dlarrayref is a number, but when no barcode file was made, dlarrayref is the actual name of the sample
	    #Ask Sacha, system command to rename a file from perl
	    qx("$mv" reads"$dlarrayref".fasta  transposon_"$dlarrayref".fasta);
	}
    }
    unlink("transposon_split_unmatched.fasta");
}#END trim_and_split

######################################################################################
# remove_transposon
# This subroutine will remove the transposon site from the splitted files. Input is 
# named transposon_$sample.fasta, output must be named $sample.fasta

sub remove_transposon {
    
    my @samples=@_;
    my $sample;
    
    foreach $sample (@samples) {
	stat("long_transposon_".$sample.".fasta");
        if (-e _) {
    	    print "Long transposon sequence given for file ".$sample.".fasta. This sequence will be trimmed using pass. \n";
    	    qx("$pass" -i transposon_"$sample".fasta -d long_transposon_"$sample".fasta -end 1 -flc 0 -fid 20 -cpu 4 2>/dev/null |"$grep" -v "^\$" > "$sample".fasta );
    	    unlink("transposon_".$sample.".fasta");
        }
        else {
    	    stat("short_transposon_".$sample.".txt");	
    	    if (-e _) {	
    		print "Short transposon sequence given for file ".$sample.".fasta. This sequence will be trimmed using fastx_barcode_splitter_trimmer. \n";
    		qx("$cat" transposon_"$sample".fasta |$fastx_split --bcfile short_transposon_"$sample".txt --prefix "" --suffix .fasta --eol --mismatches 3 2>>remove_transposon_fastxsplit.error);
    		unlink("transposon_".$sample.".fasta");
    		unlink("unmatched.fasta");
    	    }
    	    else {
    		qx("$mv" transposon_"$sample".fasta "$sample".fasta 2>>remove_transposon_mv.error );
    	    }
        } 
    }
}

# Sacha
# my @arr ;
# $arr_ref = doiets(\@arr) ;
# @arr = @{ $arr_ref } ;
#

# END remove_transposon


#######################################################################################
# pass_alignments
# This subroutine executes the pass commands to align the reads to the ta sites. Replace 
# this with your favorite aligner subroutine. input is $sample.fasta, output of the aligner must be names $tacounts

sub pass_alignments {
    my @samples=@_;
    my $sample;
    foreach $sample (@samples) {
	stat ("".$sample.".fasta");
	unless (-e _) {general::error("sample ".$sample.".fasta was not created. Your file was not split on your barcode or there was no match with your transposon sequence. Try checking your barcode and run the software without transposon sequence \n");}
	
	my $alltacounts = "allta_".$sample.".counts.txt";
	my $tacounts = "ta_".$sample.".counts.txt";
	my $wiggle = "ta_".$sample.".wiggle";
	
	#create header for ta counts file
        open (TACOUNTS, "> $tacounts");
        print TACOUNTS "counts\tnames\n";
        close (TACOUNTS);
	                	        
	#align reads of samples to the ta sites. use of cut,sort,uniq,sed and tr is a bit of a hack
	print "Aligning ".$sample.".fasta to the insertion sites using pass\n"; 
	qx("$pass" -i "$sample".fasta -d genome.ta -fle "$minhit" -b -uniq -flc 0 -cpu 4 -s -gff -S "$strand" 2>/dev/null |"$cut" -f 1 |"$sort" |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > $alltacounts);
	qx("$cat" "$alltacounts" |"$sort" -n -r |head -n "$librarysize" > "$tacounts");
	unlink("".$sample.".fasta");
		
	#awful hack. plan to incorporate the removal of the - sign and the sorting in the perl script below. GNU sort is quicker though..
	qx("$cat" "$tacounts" |"$grep" -v names|"$tr" -d '-'|"$sort" -n -k2 >$wiggle);          
	print "Counting read frequencies of genes and TA sites of ".$sample.".fasta\n";
	# rapidly add all ta site counts (the wiggle track of the genome) per gene, which positions are annotated in genome.ptt generated by gbk2fnaptt.
	my $featureend;
	my $featurestart;
	my $featureline;
	my $wiggleline;
	my $wigglepos;
	my $featurecount;
	my $featurename;
	my $wiggledepth;
	my @splet;
	my @splyt;
	my $wigglepos_remember=0;
	my $wiggledepth_remember=0;
	my $countsgene = "gene_".$sample.".counts.txt";
    
	
	open (COUNTSGENE, "> $countsgene" ) or general::error("Can't open file: ".$countsgene."");
	print COUNTSGENE "counts\tnames\n";
	open(LOCATIONS, 'genome.ptt') or general::error("Can't open file: genome.ptt");
	open(WIGGLE, "< $wiggle" ) or general::error("Can't open file: ".$wiggle."");
	while (<LOCATIONS>){
    	    $featureline = $_ ;
    	    chomp $featureline ;
    	    @splyt = split '\t', $featureline;
    	    $featurename =  $splyt[0];
    	    $featurestart = $splyt[1];
    	    $featureend = $splyt[2];
    	    $featurecount=0;
    	    next if ($featurename eq "Locus");
    	    while (<WIGGLE>){
        	$wiggleline = $_ ;
        	chomp $wiggleline ;
        	@splet = split '\t', $wiggleline;
        	$wigglepos = $splet[1];
        	$wiggledepth = $splet[0];
        	
        	# remember the position and depth for the next line in LOCATIONS
        	if ($wigglepos > $featureend) {
        	    $wigglepos_remember = $wigglepos;
        	    $wiggledepth_remember = $wiggledepth;
        	}
        	
        	#IMPORTANT! if the previous wigglepos was inside the current gene, add these counts as well. Else the counts of the first TA site will not get added to $featurecount 
        	if ($wigglepos_remember >= $featurestart && $wigglepos_remember <= $featureend) {
        	    $featurecount=$featurecount+$wiggledepth_remember;
        	    $wiggledepth_remember = 0;
        	}
        	
        	#skip all next statements. This prints the featurecount to COUNTSGENE
        	last if ($wigglepos > $featureend);
        	  
        	if ($wigglepos >= $featurestart && $wigglepos <= $featureend) {
        	    $featurecount=$featurecount+$wiggledepth;
        	    $wiggledepth_remember = 0;
        	}
    	    }
    	    print COUNTSGENE "$featurecount\t$featurename\n";
	}
	close(WIGGLE);
	close(LOCATIONS);
	close(COUNTSGENE);
    }
}# End pass_alignments


#######################################################################################
# doTAcounts
# This rubroutine prepares unique TA site counts per gene for essential gene statistics

sub doTAcounts {

    my $featureend;
    my $featurestart;
    my $featureline;
    my $wiggleline;
    my $wigglepos;
    my $featurecount;
    my $featurename;
    my $wiggledepth;
    my @splet;
    my @splyt;
    my $countsgene = "tavsgenes.txt";
    my $wiggle = "tavsta.wiggle";    
    
    print "Counting number of TA sites per gene\n";
    #perform pass alignment command, sort command and uniq command, do some magic with sed and tr to make it tab delimited.  
    qx("$pass" -i genome.ta -d genome.ta -fle 14 -b -uniq -flc 0 -cpu 4 -s -gff -S 0 2>/dev/null |"$cut" -f 1 |"$sort" |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > tavsta.txt);
    #awful hack. plan to incorporate the removal of the - sign in the perl script below. I need to sort the file anyway.
    qx("$cat" tavsta.txt |"$grep" -v names|"$tr" -d '-'|"$sort" -n -k2 >tavsta.wiggle);          
	
    # rapidly add all ta site counts (the wiggle track of the genome) per gene, which positions are annotated in genome.ptt.
    
    open (COUNTSGENE, "> $countsgene" ) or general::error("Can't open file: ".$countsgene."");
    print COUNTSGENE "counts\tnames\n";
    open(LOCATIONS, 'genome.ptt') or general::error("Can't open file: genome.ptt");
    open(WIGGLE, "< $wiggle" )  or general::error("Can't open file: ".$wiggle."");
    while (<LOCATIONS>){
        $featureline = $_ ;
        chomp $featureline ;
        @splyt = split '\t', $featureline;
        $featurename =  $splyt[0];
        $featurestart = $splyt[1];
        $featureend = $splyt[2];
        $featurecount=0;
        next if ($featurename eq "Locus");
        while (<WIGGLE>){
            $wiggleline = $_ ;
            chomp $wiggleline ;
            @splet = split '\t', $wiggleline;
            $wigglepos = $splet[1];
            $wiggledepth = $splet[0];
            last if ($wigglepos > $featureend);
            if ($wigglepos >= $featurestart && $wigglepos <= $featureend) {
            $featurecount=$featurecount+$wiggledepth;
            }
        }
        print COUNTSGENE "$featurecount\t$featurename\n";
    }
    close(WIGGLE);
    close(LOCATIONS);
    close(COUNTSGENE);

    open (TARGETSESSENTIALGENES, '>>essentialgenes_targets.txt');
    print TARGETSESSENTIALGENES "tavsgenes.txt\tcontrol\tcontrol\tinsertionsites\n";
    close (TARGETSESSENTIALGENES);
}# END doTAcounts

#######################################################################################
# doTAannotation
# This rubroutine annotates the TA sites

sub doTAannotation {

    my $featureend;
    my $featurestart;
    my $featureline;
    my $wiggleline;
    my $wigglepos;
    my $featurecount;
    my $featurename;
    my $wiggledepth;
    my $feature3;
    my $feature4;
    my $feature5;
    my $feature6;
    my $feature7;
    my $wigglename;
    my @splet;
    my @splyt;
    my $annotation = "ta.ptt";
    my $wiggle = "tavsta.wiggle3";    

    print "Annotating TA sites\n";
    #blegh. ugly hack. fix asap
    qx("$cat" tavsta.txt |"$tr" -d '-' >tavsta.wiggle2);
    qx("$paste" tavsta.wiggle2 tavsta.txt |"$sort" -n -k 2 >tavsta.wiggle3);
    # rapidly annotate the TA sites (the wiggle track of the genome) with genes, which positions are annotated in genome.ptt.
    
    open (ANNOTATION, "> $annotation" )	or general::error("Can't open file: ".$annotation."");
    print ANNOTATION "counts\tnames\n";
    open(LOCATIONS, 'genome.ptt') 	or general::error("Can't open file: genome.ptt");
    open(WIGGLE, "< $wiggle" ) 		or general::error("Can't open file: ".$wiggle."");
    while (<LOCATIONS>){
        $featureline = $_ ;
        chomp $featureline ;
        @splyt = split '\t', $featureline;
        $featurename =  $splyt[0];
        $featurestart = $splyt[1];
        $featureend = $splyt[2];
        $feature3 = $splyt[3]; 
        $feature4 = $splyt[4];
        $feature5 = $splyt[5];
        $feature6 = $splyt[6];
        $feature7 = $splyt[7];
        $featurecount=0;
        next if ($featurename eq "Locus");
        while (<WIGGLE>){
            $wiggleline = $_ ;
            chomp $wiggleline ;
            @splet = split '\t', $wiggleline;
            $wigglepos = $splet[1];
            $wiggledepth = $splet[0];
            $wigglename = $splet[3];
            last if ($wigglepos > $featureend);
            if ($wigglepos >= $featurestart && $wigglepos <= $featureend) {
            print ANNOTATION "$wigglename\t$featurename\t$featurestart\t$featureend\t$feature3\t$feature4\t$feature5\t$feature6\t$feature7\n";
            $featurecount=$featurecount+$wiggledepth;
            }
        }
    }
    close(WIGGLE);
    close(LOCATIONS);
    close(ANNOTATION);
} #END doTAannotation

#######################################################################################
# doStats
# This subroutine will perform the statistics. Has to be moved to wrapper script

sub doStats {
    print "Performing normalization using loess and TMM (edgeR)\n";
    print "Performing statistical tests (edgeR) on conditionally essential genes\n";
    qx("$Rnew" --vanilla --args $normalization < "$gene_edger" 2>>R.error);
    print "Performing statistical tests (edgeR) on essential genes\n";
    qx("$Rnew" --vanilla --args $normalization < "$esgene_edger" 2>>R.error);
    print "Performing statistical tests (edgeR) on insertion sites\n";
    qx("$Rnew" --vanilla --args $normalization < "$ta_edger" 2>>R.error);

}#END doStats

######################################################################################
# generatehtml
# This subroutine htmlifies some of the output

sub generatehtml {
    
    my $output;    
    open (OUTPUT, '> output.html' ) or general::error("Can't open file: output.html");
    print OUTPUT <<END;
<br>
<b>Conditionally essential genes:</b><br>
<a href="$sessiondir/gene_MDS_non_normalized.png" target="_blank" >MDS plot of non\* normalized read counts on genes</a><br>
<a href="$sessiondir/gene_MDS_normalized.png " target="_blank" >MDS plot of normalized read counts on genes</a><br>
<a href="$sessiondir/gene_plot.png" target="_blank" >Plot of signal versus log2 fold change on genes</a><br>
<a href="$sessiondir/gene_alloutputmerged.tsv" target="_blank" >Combined table of conditionally essential genes (.tsv)</a><br><br>
Add conditionally essential gene data to database for MINOMICS visualization.<br>
END
    $output=0;
    $output = add2db_link::add2db_html_link('madatabase','gene_output_add2db.html',$sessiondir.'gene_output_add2db.html',$workdir.'gene_minomics.tsv','madatabase','Multiple','1');
    print OUTPUT " ".$output." \n";

    print OUTPUT <<END;
<hr>
<b>Conditionally essential insertion sites:</b><br>
<a href="$sessiondir/ta_MDS_non_normalized.png" target="_blank" >MDS plot of non\* normalized read counts on insertion sites</a><br>
<a href="$sessiondir/ta_MDS_normalized.png" target="_blank" >MDS plot of normalized read counts on insertion sites</a><br>
<a href="$sessiondir/ta_plot.png" target="_blank" >Plot of signal versus log2 fold change on insertion sites</a><br>
<a href="$sessiondir/ta_alloutputmerged.tsv" target="_blank" >Combined table of conditionally essential insertion sites (.tsv)</a><br><br>
END
    print OUTPUT <<END;
<hr>
<b>Essential genes.</b> <br>
<a href="$sessiondir/essentialgenes_MDS_non_normalized.png " target="_blank" >MDS plot of non\* normalized read counts on essential genes</a><br>
<a href="$sessiondir/essentialgenes_MDS_normalized.png" target="_blank" >MDS plot of normalized data read counts on essential genes</a><br>
<a href="$sessiondir/essentialgenes_plot.png " target="_blank" >Plot of signal versus log2 fold change on essential genes</a><br>
<a href="$sessiondir/essentialgenes_alloutputmerged.tsv" target="_blank" >Combined table for essential genes (.tsv)</a><br><br>
Add essential gene data to database for MINOMICS visualization.<br>
END
    $output = 0;
    $output = add2db_link::add2db_html_link('madatabase','essentialgenes_output_add2db.html',$sessiondir.'essentialgenes_output_add2db.html',$workdir.'essentialgenes_minomics.tsv','madatabase','Multiple','1');
    print OUTPUT " ".$output." \n";

    print OUTPUT <<END;
<hr>
<b>Essential genes with genomic location bias removed using loess normalization</b><br>
<a href="$sessiondir/loess_essentialgenes_MDS_non_normalized.png" target="_blank" >MDS plot of non\* normalized read counts on essential genes</a><br>
<a href="$sessiondir/loess_essentialgenes_MDS_normalized.png" target="_blank" >MDS plot of normalized data read counts on genes</a><br>
<a href="$sessiondir/loess_essentialgenes_plot.png " target="_blank" >Plot of signal versus log2 fold change on essential genes</a><br>
<a href="$sessiondir/loess_essentialgenes_alloutputmerged.tsv" target="_blank" >Combined table for essential genes (.tsv)</a><br><br>
Add essential gene data with genomic location bias removed to database for MINOMICS visualization.<br>
END
    $output=0;
    $output = add2db_link::add2db_html_link('madatabase','loess_essentialgenes_output_add2db.html',$sessiondir.'loess_essentialgenes_output_add2db.html',$workdir.'loess_essentialgenes_minomics..tsv','madatabase','Multiple','1');
    print OUTPUT " ".$output." \n";
    
    print OUTPUT <<END;
<hr>
<a href="http://bamics2.cmbi.ru.nl/websoftware/minomics" target="_blank" > Run MINOMICS on database data for visualization. Remember to add the data first. </a>
<br>
<br>
<hr>
\* Correction based on library size always occurs. This is inherent to the statistical test used. <br>


END
    close OUTPUT ;    
}#END generatehtml
