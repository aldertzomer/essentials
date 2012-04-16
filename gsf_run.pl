#!/usr/local/bin/perl
# GSF perl script. Calls scripts for parsing the genbank (module 1)  and generating the insertion site flanking sequences (module 2), 
# parses configfile, downloads reads, does alignments  (module 3)  and runs statistical analysis module 4 (R)
# prerequisites: R, GNU tools, curl, pass. fq_all2std_noqual.pl, fastx_barcode_splitter_trimmer, gbk2pttfna.pl, tasite.pl
# fastx_barcode_splitter_trimmer requires Text::LevenshteinXS if you want to split your barcodes quickly

use strict;
use warnings;
use fgweb::general::general;
#use fgweb::general::add2db_link;
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
my $barcode_end = $ENV{'BARCODE_END'};
my $transposon_end = $ENV{'TRANSPOSON_END'};
my $minhit = $ENV{'MINHIT'};
my $tasitelength = $ENV{'TASITELENGTH'};
my $strand = $ENV{'STRAND'};
my $librarysize = 2 * $ENV{'LIBRARYSIZE'};
my $statmethod = $ENV{'STATMETHOD'};
my $adjustment = $ENV{'ADJUSTMENT'};
my $dispersion = $ENV{'DISPERSION'};
my $smoothing = $ENV{'SMOOTHING'};
my $ptt = $ENV{'GENETRUNCATION'};
my $loess = $ENV{'LOESS'};
my $repeat = $ENV{'REPEAT'};
my $uniqoption ;
my $minomicsreads = $ENV{'MINOMICSREADS'};

my $zip = $ENV{'DOZIP'};


if ($repeat eq "No") {$uniqoption = " ";}
if ($repeat eq "Yes") {$uniqoption = "-uniq";}

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

#fix:
my $revseq = '/usr/local/bin/revseq';

#perl scripts for parsing the genbank, the insertion site finder, the genome splitter, the fastx tool kit and maq
#These script are either written by someone else and modified, or based on existing code.

my $gbk2fna 		= $vars->{cgi}{root}.'gsf/gbk2pttfna.pl';
my $gbk2fna3             = $vars->{cgi}{root}.'gsf/gbk2pttfna3.pl';
my $tasitescript 	= $vars->{cgi}{root}.'gsf/tasite.pl';
my $split		= $vars->{cgi}{root}.'gsf/split.pl';
my $fq_all2std  	= $vars->{cgi}{root}.'gsf/fq_all2std_noqual.pl';
my $fastx_split 	= $vars->{cgi}{root}.'gsf/fastx_barcode_splitter_trimmer.pl';

# R scripts for stats

my $loessscript		= $vars->{cgi}{root}."gsf/R/loess.R";
my $loessgenescript	= $vars->{cgi}{root}."gsf/R/loess_genes.R";
my $gene_edger 		= $vars->{cgi}{root}."gsf/R/".$statmethod."_gene_edgeR.R";
my $esgene_edger 	= $vars->{cgi}{root}."gsf/R/".$statmethod."_essentialgenes_edgeR.R";
my $ta_edger 		= $vars->{cgi}{root}."gsf/R/".$statmethod."_ta_edgeR.R";

use lib "/home/fgweb/cgi-bin/gsf/";
use add2db_link;

my $archivefile = 'archive.zip';

######################################################################################
# Run programs down here. Comment out sections to speed up a rerun or to debug.

#cleanup previous run
unlink("allinsertions.txt");
unlink("samplenames.txt");
unlink (<*.error>);

# Parse FGweb generated configfile if exists.
my $fgwebconfigfile = "configfile_uploaded.txt";
stat($configfile);
unless (-e _) {
    print "Generating downloadable configfile\n";
    my $fgweblineconfig;
    my $fgwebvarname;
    my $fgwebvarcontent;
    my @fgwebsplot;

    open(FGWEBCONFIGFILE, $fgwebconfigfile) or die ("Can't open file:" .$fgwebconfigfile);
    open (CONFIGFILE, "> $configfile");
    print CONFIGFILE "Link\tbarcode\ttransposon_end\tsampletype\tlibrary\tsampleformat\tcompression\tname\n";
    
    while (<FGWEBCONFIGFILE>){
	$fgweblineconfig = $_ ;
	chomp $fgweblineconfig ;
        @fgwebsplot = split '\=', $fgweblineconfig;
        $fgwebvarname = $fgwebsplot[0] ;
        $fgwebvarcontent = $fgwebsplot[1] ;
        next if ($fgwebvarname eq "config_experimentname");
        next if ($fgwebvarname eq "");
        next if ($fgwebvarname eq " ");
        if ($fgwebvarname=~ m/^config_name/) {print CONFIGFILE "$fgwebvarcontent\n";}
        else {print CONFIGFILE "$fgwebvarcontent\t";}
	}
    close FGWEBCONFIGFILE;
    close CONFIGFILE;
}


# Parsing genbank
print "Parsing selected genbank file. \n";
qx("$perl" "$gbk2fna3" -gbk="$genbank" -s=genome.fasta -sa=genome_all.ptt 2>> gbk2fna.error);

# sorting ptt file with genes and proteins on genename and then be annotation length
my @pttdata;
my @ordered_pttdata;
open(ALLPTT, "genome_all.ptt") or general::error("Can't open file: genome_all.ptt");
open (ALLPTT2, "> genome_all_sorted.ptt");

while (<ALLPTT>){
    chomp;
    push(@pttdata, $_);
} 

@ordered_pttdata = sort { (split '\t', $a)[0]  cmp (split '\t', $b)[0] || (length $b  <=> length $a) } @pttdata;
print ALLPTT2 join( "\n", @ordered_pttdata )."\n";
close ALLPTT;
close ALLPTT2;


# TODO: get rid of crap below and do the duplicate search on the array.
my $alllocus;
my $allstart;
my $allstop ;
my $allstrand ;
my $alllength ;
my $allpid ;
my $allgene ;
my $allproduct ;
my $alllocusold ;
my $allpttline;
my @splitallptt;
$alllocusold = 0;

open(ALLPTT3, "genome_all_sorted.ptt") or general::error("Can't open file: genome_all_sorted.ptt");
open (PTT, "> genome.ptt");

while (<ALLPTT3>){
    $allpttline = $_ ;
    chomp $allpttline ;
    @splitallptt = split '\t', $allpttline;
    $alllocus = $splitallptt[0] ;
    $allstart = $splitallptt[1] ;
    $allstop = $splitallptt[2] ;
    $allstrand = $splitallptt[3] ;
    $alllength =  $splitallptt[4] ;
    $allpid = $splitallptt[5];
    $allgene = $splitallptt[6];
    $allproduct = $splitallptt[7];
    if ($alllocus ne $alllocusold) {
        print PTT "$alllocus\t$allstart\t$allstop\t$allstrand\t$alllength\t$allpid\t$allgene\t$allproduct\n";
    }
    $alllocusold = $alllocus;
}
close PTT;
close ALLPTT3;

# Truncate end of genes to filter out function-retaining C-terminal transposon insertions
truncategenes();

# Finding TA insertion site flanking sequences. 
if ($insertion eq "TA") {
    print "Searching TA sites and indexing insertion site flanking sequences ($tasitelength bp) of selected genome\n";
    qx("$perl" "$tasitescript" genome.fasta "$insertion" "$tasitelength" >genome.ta 2>> tasitescript.error);
}

# Fragmenting genomic DNA for detection of non-unique insertion site flanking sequences.
if ($insertion eq "random") {
    print "Indexing insertion site flanking sequences ($tasitelength bp) of selected genome\n";
    qx("$perl" "$split" genome.fasta "$tasitelength" >genome.ta_fw 2>>splitscript.error);
    qx("$revseq" -sequence genome.ta_fw -outseq genome.ta_rev);
    qx("$cat" genome.ta_fw genome.ta_rev > genome.ta);
}

# Run the analysis
get_and_convert();

# Run the statistics
doStats();

# Generate html
generatehtml();

# zip everything 
# bamics2 still has an old perl version. cannot use a bare readdir in a while loop and set $_ on every iteration. have to use @

if ($zip eq "Yes") {
    
    use Archive::Zip qw( :ERROR_CODES :CONSTANTS );
    
    my $archive = Archive::Zip->new();
    my $dirfile ;
    
    unlink ($archivefile);
    opendir DIR, $workdir or die "cannot open dir $workdir: $!";
    my @files= readdir DIR;
    closedir DIR;
    foreach $dirfile (@files) {if ($dirfile ne "."){if ($dirfile ne ".."){$archive->addFile($dirfile);}}}
    unless ( $archive->writeToFileNamed($archivefile) == AZ_OK ) {     
        die 'write error';
    }

}

# debug option
# qx(chmod 777 -R "$workdir/");

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
    my $samplename;
    my $decompress;
    my $library ;
    my %config ;
    my @splot ;
    my @splut ;
    #unused    my $transposonmismatch;

    #print "Preparing headers of output files\n";
    writeheaders();
    
    # parse the configfile line by line. each line contains the link, the barcode, the transposon sequence, the sampletype (control or target)  and the input format
    # If the barcode or transposon is not present, use N as sequence. 
    
    print "\nParsing configfile\n";
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
        $samplename =  $splot[7];
        
	next if ($inputlink eq "Link");
        next if ($inputlink eq "link");
        next if ($inputlink eq "#");
        next if ($inputlink eq "#link");
        next if ($inputlink eq "#Link");
	next if ($inputlink eq "");
	next if ($inputlink eq " ");
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
	    if ($inputformat eq "fasta") { qx("$curl" -f --compressed "$inputlink" --stderr "$readfile".error $decompress > "$readfile");}
	    if ($inputformat eq "fastq") { qx("$curl" -f --compressed "$inputlink" --stderr "$readfile".error $decompress |"$fq_all2std" fq2fa   	> "$readfile");}
	    if ($inputformat eq "scarf") { qx("$curl" -f --compressed "$inputlink" --stderr "$readfile".error $decompress |"$fq_all2std" scarf2std	|"$fq_all2std" fq2fa  > "$readfile");}
	    if ($inputformat eq "export"){ qx("$curl" -f --compressed "$inputlink" --stderr "$readfile".error $decompress |"$fq_all2std" export2std	|"$fq_all2std" fq2fa  > "$readfile");}
	    if ($inputformat eq "solid") { qx("$curl" -f --compressed "$inputlink" --stderr "$readfile".error $decompress |"$fq_all2std" csfa2std	|"$fq_all2std" fq2fa  > "$readfile");}

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
            
            #create samplename file for later use
	    if ($sampletype ne "ignore") {
		open (SAMPLENAME, ">> samplenames.txt");
		if ($samplename ne "") {print SAMPLENAME "".$samplename."_".$samplecount."_".$sampletype."\n";}
		if ($samplename eq "") {print SAMPLENAME "sample".$samplecount."_".$sampletype."\n";}
        	close (SAMPLENAME);
            }
            
            #if there's no barcode, there's also no need to split. Immediately give the readfile a sample name instead of a number                                  
	    if ($barcode eq "N") {
		    print "File $readfile has no barcode assigned. Assuming the entire file is sample $sample\n";
		    countlines($readfile);
		    qx("$mv" "$readfile" transposon_"$sample".fasta 2>> get_and_convert_mv.error);
	    }	    
	    
	    #create transposon site file for later use (long or short, has an effect on downstream processing)
	    if (length $transposon < 13) {
    		if ($transposon ne "N") {
    		    my $transposonfile = "short_transposon_".$sample.".txt";
		    open (TRANSPOSON, "> $transposonfile");
		    print TRANSPOSON "$sample\t$transposon\n";
		    close (TRANSPOSON);
		}
	    }
	    if (length $transposon >= 13) {
    		my $transposonfile = "long_transposon_".$sample.".fasta";
		open (TRANSPOSON, "> $transposonfile");
		print TRANSPOSON "\>transposon\n";
		print TRANSPOSON "$transposon\n";
		close (TRANSPOSON);
	    }
	    
	    #add sample to file with control and target samples for gene edgeR
	    if ($sampletype ne "ignore") {
		my $countsgene = "gene_".$sample.".counts.txt";
		open (TARGETSGENES, ">> gene_targets.txt");
		print TARGETSGENES "".$countsgene."\t".$sampletype."\t".$sampletype."\t".$library."\n";
		close (TARGETSGENES);
	    
		#add sample to file with all samples as control samples for essential genes edgeR
		open (TARGETSESSENTIALGENES, ">> essentialgenes_targets.txt");
		print TARGETSESSENTIALGENES "".$countsgene."\ttarget\ttarget\t".$library."\n";
		close (TARGETSESSENTIALGENES);
		
		#add sample to file with control and target samples for insertion site edgeR
		my $tacounts = "ta_".$sample.".counts.txt";
		open (TARGETSTA, ">> ta_targets.txt");
		print TARGETSTA "".$tacounts."\t".$sampletype."\t".$sampletype."\t".$library."\n";
		close (TARGETSTA);
		push(@samples, $sample);
	    }
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
	    
            #Add samplename to file for later use
            if ($sampletype ne "ignore") {
                open (SAMPLENAME, ">> samplenames.txt");
                if ($samplename ne "") {print SAMPLENAME "".$samplename."_".$samplecount."_".$sampletype."\n";}
                if ($samplename eq "") {print SAMPLENAME "sample".$samplecount."_".$sampletype."\n";}
                close (SAMPLENAME);
            }

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
	    if ($sampletype ne "ignore") {
		my $countsgene = "gene_".$sample.".counts.txt";
		open (TARGETSGENES, ">> gene_targets.txt");
		print TARGETSGENES "".$countsgene."\t".$sampletype."\t".$sampletype."\t".$library."\n";
		close (TARGETSGENES);
	    
		#add sample to file with all samples as control samples for essential genes edgeR
		open (TARGETSESSENTIALGENES, ">> essentialgenes_targets.txt");
		print TARGETSESSENTIALGENES "".$countsgene."\ttarget\ttarget\t".$library."\n";
		close (TARGETSESSENTIALGENES);
		
		#add sample to file with control and target samples for insertion site edgeR
		my $tacounts = "ta_".$sample.".counts.txt";
		open (TARGETSTA, ">> ta_targets.txt");
		print TARGETSTA "".$tacounts."\t".$sampletype."\t".$sampletype."\t".$library."\n";
		close (TARGETSTA);
		push(@samples, $sample);
	    }
	}
    }
    #enter the numbers of the downloaded files into the trimming and splitting routine creating per sample a fasta file
    print "\n";
    trim_and_split(@downloadedfiles);
    print "\n";
    
    #remove the long or short transposon if necessary
    remove_transposon(@samples);
    print "\n";
    
    # align the fasta files of the samples to the insertion site flanking sequences, count the number of reads per insertion site and count the number of reads per gene
    pass_alignments(@samples);

    #calculate number of unique insertion site flanking sequences per gene
    doTAcounts();
    
    #annotate the insertion site flanking sequences
    doTAannotation();
    print "\nCounting total number of insertion site flanking sequences hit\n";
    countlines2("allinsertionsuniq.txt");
        
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
	countlines("reads".$dlarrayref.".fasta");
    	print "Splitting reads".$dlarrayref.".fasta\n";
    	stat("barcodes".$dlarrayref.".txt");
    	if (-e _) {
	    qx("$cat" reads"$dlarrayref".fasta | $fastx_split --bcfile barcodes"$dlarrayref".txt --prefix transposon_split_ --suffix .fasta --"$barcode_end" --mismatches "$barcodemismatch"  2>>trim_and_split_fastxsplit.error);
	    unlink("reads".$dlarrayref.".fasta");
	    unlink (<*ignore*>); #deletes all samples that are not used. 
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
    	    print "Long transposon sequence given for file ".$sample.".fasta. This sequence will be filtered and trimmed using pass. \n";
    	    print "Counting number of sequences before trimming\n";
	    countlines("transposon_".$sample.".fasta");
	    if ($transposon_end eq "bol") { #awful hack. but pass is so much faster at removing transposon ends, and it supports variable lengths of the transposon. Unfortunately it does not support tag removal at the beginning of the line. Therefore we have to revseq
		qx("$revseq" -sequence transposon_"$sample".fasta -outseq transposon_"$sample".rev);
		unlink("transposon_".$sample.".fasta");
		qx("mv" transposon_"$sample".rev transposon_"$sample".fasta);
		#yes also the transposon end needs to be revseq'ed.
		qx("$revseq" -sequence long_transposon_"$sample".fasta -outseq long_transposon_"$sample".rev);
		unlink("long_transposon_".$sample.".fasta");
		qx("mv" long_transposon_"$sample".rev  long_transposon_"$sample".fasta);
	    }
	    qx("$pass" -i transposon_"$sample".fasta -d long_transposon_"$sample".fasta -end 1 -fid 20  -cpu 4 2>/dev/null |"$grep" -v "^\$" > "$sample".fasta );
    	    if ($transposon_end eq "bol") { #revseq it back to have the original situation again
                qx("$revseq" -sequence "$sample".fasta -outseq "$sample".rev);
                unlink("".$sample.".fasta");
                qx("mv" "$sample".rev "$sample".fasta);
		print "Counting number of sequences after trimming\n";
        	countlines("".$sample.".fasta");
	    }
	    unlink("transposon_".$sample.".fasta");
        }
        else {
    	    stat("short_transposon_".$sample.".txt");	
    	    if (-e _) {	
    		print "Counting number of sequences before trimming\n";
    		countlines("transposon_".$sample.".fasta");
    		print "Short transposon sequence given for file ".$sample.".fasta. This sequence will be trimmed using fastx_barcode_splitter_trimmer. \n";
    		qx("$cat" transposon_"$sample".fasta |"$fastx_split" --bcfile short_transposon_"$sample".txt --prefix "" --suffix .fasta --"$transposon_end" --mismatches 2 2>>remove_transposon_fastxsplit.error);
    		unlink("transposon_".$sample.".fasta");
    		unlink("unmatched.fasta");
    		print "Counting number of sequences after trimming\n";
    		countlines("".$sample.".fasta");
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
# This subroutine executes the pass commands to align the reads to the insertion site flanking sequences or to the genome.  
# input is $sample.fasta, output of the aligner is named $tacounts

sub pass_alignments {
    my @samples=@_;
    my $sample;
    foreach $sample (@samples) {
	stat ("".$sample.".fasta");
	unless (-e _) {general::error("sample ".$sample.".fasta was not created. Your file was not split on your barcode or there was no match with your transposon sequence. Try checking your barcode and run the software without transposon sequence \n");}
	
	my $alltacounts = "allta_".$sample.".counts.txt";
	my $tacounts = "ta_".$sample.".counts.txt";
	my $passout = "passout_".$sample.".txt";
	my $wiggle = "ta_".$sample.".wiggle";
	
	#create header for ta counts file
        open (TACOUNTS, "> $tacounts");
        print TACOUNTS "counts\tnames\n";
        close (TACOUNTS);

	if ($insertion eq "TA") {                	        
	    #align reads of samples to the ta sites. use of cut,sort,uniq,sed and tr is a bit of a hack
	    print "\nAligning ".$sample.".fasta to the insertion site flanking sequences using pass\n";
	    countlines("".$sample.".fasta");
	    print "pass -i ".$sample.".fasta -d genome.ta -fle ".$minhit." -b ".$uniqoption."  -cpu 4 -j -g 3 -s -gff -S ".$strand." \n"; 
	    qx("$pass" -i "$sample".fasta -d genome.ta -fle "$minhit" -b "$uniqoption"   -cpu 4 -j -g 3 -s -gff -S "$strand" 2>/dev/null |"$cut" -f 1 |"$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > $alltacounts);
	    unlink("".$sample.".fasta");
	}
	if ($insertion eq "random") {
	    #align reads of samples to the genome. use of cut,sort,uniq,sed and tr is a bit of a hack
	    print "\nAligning ".$sample.".fasta to the genome using pass\n"; 
	    countlines("".$sample.".fasta");
	    print "".$pass." -i ".$sample.".fasta -d genome.fasta -fle ".$minhit." -b ".$uniqoption."  -j -g 3 -cpu 4 -s -gff -S 2 \n";
	    qx("$pass" -i "$sample".fasta -d genome.fasta -fle "$minhit" -b "$uniqoption"  -j -g 3 -cpu 4 -s -gff -S 2  2>/dev/null > $passout);
	    if ($transposon_end eq "bol") {
		qx("$cat" "$passout" |"$grep" "	-	" |"$cut" -f 5 | "$sed" 's?^?-?' |"$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > $alltacounts);
		qx("$cat" "$passout" |"$grep" "	+	" |"$cut" -f 4 | "$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' >> $alltacounts);
		}
	    if ($transposon_end eq "eol") {
		qx("$cat" "$passout" |"$grep" "	+	" |"$cut" -f 5 | "$sed" 's?^?-?' |"$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > $alltacounts);
		qx("$cat" "$passout" |"$grep" "	-	" |"$cut" -f 4 | "$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' >> $alltacounts);
		}    
	    
	    unlink("".$sample.".fasta");
	    unlink("".$passout."");
	}
	
	countwiggle($alltacounts);
	print "Retaining only the $librarysize insertion site flanking sequences  with the highest numbers of reads\n";
	
	qx("$cat" "$alltacounts" |"$sort" -n -r |head -n "$librarysize" > "$tacounts");
    	
    	if ($loess eq "Yes") {
	    print "Loess normalizing ".$tacounts." for read count bias\n";
	    qx("$Rnew" --vanilla --args "$tacounts" < "$loessscript" 2>>R.error);
	}    	         
    	                                                         
    	
    	#awful hack. plan to incorporate the removal of the - sign and the sorting in the perl script below. GNU sort is quicker though..
	qx("$cat" "$tacounts" >>allinsertions.txt);
	qx("$cat" "$tacounts" |"$grep" -v names|"$tr" -d '-'|"$sort" -n -k2 >$wiggle);          
	qx("$cat" $ptt |"$sort" -k2 -n > genome.ptt.tmp);
	qx("mv" genome.ptt.tmp $ptt);
	print "Counting read frequencies of genes and insertion site flanking sequences of ".$sample.".fasta\n";
	# rapidly add all insertion site counts (the wiggle track of the genome) per gene, which positions are annotated in genome.ptt generated by gbk2fnaptt.
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
	my $countsgene = "gene_".$sample.".counts.txt";
	my $lengthline ;
	my @buffer = (0) x 10000;
	my $filepos = 0 ;
	$featurecount=0;
	
	open (COUNTSGENE, "> $countsgene" ) or general::error("Can't open file: ".$countsgene."");
	print COUNTSGENE "counts\tnames\n";
	open(LOCATIONS, $ptt) or general::error("Can't open file: ".$ptt."");
	open(WIGGLE, "< $wiggle" ) or general::error("Can't open file: ".$wiggle."");
	while (<LOCATIONS>){
    	    $featureline = $_ ;
    	    chomp $featureline ;
    	    @splyt = split '\t', $featureline;
    	    $featurename =  $splyt[0];
    	    $featurestart = $splyt[1];
    	    $featureend = $splyt[2];
    	    next if ($featurename eq "!Locus");
    	    $featurecount = 0;
    	    $wiggledepth = 0;
    	    while (<WIGGLE>){
        	$filepos = tell(WIGGLE);
        	shift @buffer;
        	push  @buffer, $filepos;
        	
        	$wiggleline = $_ ;
        	#$lengthline = $_;
        	chomp $wiggleline ;
        	@splet = split '\t', $wiggleline;
        	$wigglepos = $splet[1];
        	$wiggledepth = $splet[0];
    		if ($featurestart <= $wigglepos && $featureend >= $wigglepos) {$featurecount = $featurecount + $wiggledepth;}
    		if ($featureend < $wigglepos) {last}
    	    }
    	    print COUNTSGENE "$featurecount\t$featurename\n";
    	    seek(WIGGLE, $buffer[0], 0);
    	    #seek(WIGGLE, 0, 0);
    	    #alternative to speed up. Insertions sites within overlapping genes will not be counted correctly when using this
    	    #seek(WIGGLE, -length($lengthline), 1);
	}
	close(WIGGLE);
	close(LOCATIONS);
	close(COUNTSGENE);
	if ($loess eq "Yes") {
	    print "Loess normalizing ".$countsgene." for insertion frequency bias\n";
	    qx("$Rnew" --vanilla --args $countsgene < "$loessgenescript" 2>>R.error);
	}
    }
}# End pass_alignments


#######################################################################################
# doTAcounts
# This rubroutine prepares unique insertion site counts per gene for essential gene statistics

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
    my $wigglepos_remember=0;
    my $wiggledepth_remember=0;
                   
    print "\nCounting number of possible unique insertion site flanking sequences per gene\n";
    #perform pass alignment command, sort command and uniq command, do some magic with sed and tr to make it tab delimited.  

    if ($insertion eq "TA") { 
	print "\nCounting number of unique insertion site flanking sequences in genome.ta before removal of non-uniques \n";
	countlines("genome.ta");
	qx("$pass" -i genome.ta -d genome.ta -b "$uniqoption" -flc 0 -cpu 4 -gff -S 0 2>/dev/null |"$cut" -f 1 |"$sort" |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > tavsta.txt);
	}

    if ($insertion eq "random") {
    	print "\nCounting number of unique insertion site flanking sequences in genome.ta before removal of non-uniques \n";
	countlines("genome.ta");
        print "\nDetecting non-unique insertion sites\n";
	print "".$pass." -i genome.ta -d genome.fasta -fle ".$minhit." -b ".$uniqoption."  -cpu 4 -s -gff -S 2 \n";
        qx("$pass" -i genome.ta -d genome.fasta -fle "$minhit" -b "$uniqoption"  -cpu 4 -j -g 3 -s -gff -S 2  2>/dev/null > passout_tavsta.txt);
        if ($transposon_end eq "bol") {
            qx("$cat" passout_tavsta.txt |"$grep" "	-	" |"$cut" -f 5 | "$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > tavsta.txt);
            qx("$cat" passout_tavsta.txt |"$grep" "	+	" |"$cut" -f 4 | "$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' >> tavsta.txt);
	    }
        if ($transposon_end eq "eol") {
            qx("$cat" passout_tavsta.txt |"$grep" "	+	" |"$cut" -f 5 | "$sed" 's?^?-?' |"$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > tavsta.txt);
            qx("$cat" passout_tavsta.txt |"$grep" "	-	" |"$cut" -f 4 | "$sed" 's?^?-?' |"$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' >> tavsta.txt);
	    }
        unlink("passout_tavsta.txt");

    	#qx("$pass" -i genome.ta -d genome.fasta -b "$uniqoption" -flc 0 -cpu 4 -gff -S 2 2>/dev/null |"$cut" -f 4 |"$sort" |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > tavsta.txt);
	}

    #awful hack. plan to incorporate the removal of the - sign in the perl script below. I need to sort the file anyway.
    print "Counting number of unique insertion site flanking sequences in genome.ta after removal of non-uniques \n";
    countlines2("tavsta.txt");
    qx("$cat" tavsta.txt |"$grep" -v names|"$tr" -d '-'|"$sort" -n -k2 >tavsta.wiggle);          
    qx("$cat" $ptt |"$sort" -k2 -n > genome.ptt.tmp);
    qx("mv" genome.ptt.tmp $ptt);
	
    # rapidly add all ta site counts (the wiggle track of the genome) per gene, which positions are annotated in genome.ptt.
    
    open (COUNTSGENE, "> $countsgene" ) or general::error("Can't open file: ".$countsgene."");
    print COUNTSGENE "counts\tnames\n";
    open(LOCATIONS, $ptt) or general::error("Can't open file: ".$ptt."");
    open(WIGGLE, "< $wiggle" )  or general::error("Can't open file: ".$wiggle."");
    
    my @buffer = (0) x 10000;
    my $filepos = 0 ;
    
    while (<LOCATIONS>){
        $featureline = $_ ;
        chomp $featureline ;
        @splyt = split '\t', $featureline;
        $featurename =  $splyt[0];
        $featurestart = $splyt[1];
        $featureend = $splyt[2];
        $featurecount=0;
        next if ($featurename eq "!Locus");
        while (<WIGGLE>){
    	    $filepos = tell(WIGGLE);
    	    shift @buffer;
    	    push  @buffer, $filepos;
    	    
    	    $wiggleline = $_ ;
    	    #$lengthline = $_ ;
    	    chomp $wiggleline ;
    	    @splet = split '\t', $wiggleline;
    	    $wigglepos = $splet[1];
    	    $wiggledepth = $splet[0];
	    if ($featurestart <= $wigglepos && $featureend >= $wigglepos) {$featurecount = $featurecount + $wiggledepth;}
	    if ($featureend < $wigglepos) {last}
	    
	}
	print COUNTSGENE "$featurecount\t$featurename\n";
	seek(WIGGLE, $buffer[0], 0);
	#seek(WIGGLE, 0, 0);
	#alternative to speed up. Insertions sites within overlapping genes will not be counted correctly when using this
	#seek(WIGGLE, -length($lengthline), 1);
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
# This rubroutine annotates the insertion site flanking sequences

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
    my $wigglepos_remember=0;
    my $wigglename_remember=0;
    print "\nAnnotating insertion site flanking sequences\n";
    my $lengthline;
    #ugly hack. fix asap
    qx("cat" allinsertions.txt |"$cut" -f 2|"$sort" -n |"$uniq" -c |"$sed" 's/^\ *//g' |"$tr" ' ' '\t' > allinsertionsuniq.txt);
    qx("$cat" allinsertionsuniq.txt |"$tr" -d '-' >tavsta.wiggle2);
    qx("$paste" tavsta.wiggle2 allinsertionsuniq.txt |"$sort" -n -k 2 >tavsta.wiggle3);
    qx("$cat" $ptt |"$sort" -k2 -n > genome.ptt.tmp);
    qx("mv" genome.ptt.tmp $ptt);

    # rapidly annotate the insertion site flanking sequences (the wiggle track of the genome) with genes, which positions are annotated in genome.ptt.
    
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
        next if ($featurename eq "!Locus");
        while (<WIGGLE>){
            $wiggleline = $_ ;
            $lengthline = $_;
            chomp $wiggleline ;
            @splet = split '\t', $wiggleline;
            $wigglepos = $splet[1];
            $wiggledepth = $splet[0];
            $wigglename = $splet[3];
    	    
    	    if ($featurestart <= $wigglepos && $featureend >= $wigglepos) {
    		print ANNOTATION "$wigglename\t$featurename\t$featurestart\t$featureend\t$feature3\t$feature4\t$feature5\t$feature6\t$feature7\n";
    	    }
    	    if ($featureend < $wigglepos) {last}
    	}
	seek(WIGGLE, -length($lengthline), 1);
    }
    close(WIGGLE);
    close(LOCATIONS);
    close(ANNOTATION);
} #END doTAannotation

#######################################################################################
# doStats
# This subroutine will perform the statistics.


sub doStats {
    print "\nPerforming statistical tests (edgeR) on conditionally essential genes\n";
    print "".$Rnew." --vanilla --args $normalization $adjustment $dispersion $smoothing < ".$gene_edger."\n" ;
    qx("$Rnew" --vanilla --args $normalization $adjustment $dispersion $smoothing $minomicsreads < "$gene_edger" 2>>R.error );
    print "Performing statistical tests (edgeR) on essential genes\n";
    print "".$Rnew." --vanilla --args $normalization $adjustment $dispersion $smoothing < ".$esgene_edger."\n";
    qx("$Rnew" --vanilla --args $normalization $adjustment $dispersion $smoothing < "$esgene_edger" 2>>R.error );
    print "Performing statistical tests (edgeR) on insertion site flanking sequences\n";
    print "".$Rnew." --vanilla --args $normalization $adjustment $dispersion $smoothing $minomicsreads < ".$ta_edger."\n"; 
    qx("$Rnew" --vanilla --args $normalization $adjustment $dispersion $smoothing $minomicsreads < "$ta_edger" 2>>R.error );

}#END doStats

######################################################################################
# generatehtml
# This subroutine htmlifies some of the output

sub generatehtml {
    
    my $output;    
    open (OUTPUT, '> output.html' ) or general::error("Can't open file: output.html");
    print OUTPUT <<END;
<br>
<a href="$sessiondir/configfile.txt" target="_blank" >Configuration file used for this run</a><br>
<a href="$sessiondir/run_output.txt" target="_blank" >Log produced by this run</a><br>
<hr>
END
    if (-e "gene_alloutputmerged.tsv") {
	print OUTPUT <<END;
<b>Conditionally essential genes</b><br>
<a href="$sessiondir/gene_PCA_non_normalized.png" target="_blank" >PCA plot of non normalized read counts on genes</a><br>
<a href="$sessiondir/gene_PCA_normalized.png " target="_blank" >PCA plot of normalized read counts on genes</a><br>
<a href="$sessiondir/gene_plot.png" target="_blank" >Plot of signal versus log2 fold change on genes</a><br>
<a href="$sessiondir/gene_densityplot.png " target="_blank" >Density Plot of log2 fold change with putative fold change cutoff(s)</a><br>
<a href="$sessiondir/gene_alloutputmerged.tsv" target="_blank" >Combined table of conditionally essential genes (.tsv)</a><br><br>
Add conditionally essential gene data to database for comparison of different experiments (also required for MINOMICS visualization).<br>
END
	$output=0;
	$output = add2db_link::add2db_html_link('madatabase','gene_output_add2db.html',$sessiondir.'gene_output_add2db.html',$workdir.'gene_minomics.tsv','madatabase','Multiple','1','resource=notappl','group1=none','group2=none','basenumber=2');
	print OUTPUT " ".$output." \n";
    } else {print OUTPUT "No conditionally essential genes analysis performed. Did you include target and control samples? <br> <hr> \n";}
    
    if (-e "ta_alloutputmerged.tsv") {
        print OUTPUT <<END;
<hr>
<b>Conditionally essential insertion sites (flanking sequences)</b><br>
<a href="$sessiondir/ta_PCA_non_normalized.png" target="_blank" >PCA plot of non normalized read counts on insertion site flanking sequences</a><br>
<a href="$sessiondir/ta_PCA_normalized.png" target="_blank" >PCA plot of normalized read counts on insertion site flanking sequences</a><br>
<a href="$sessiondir/ta_plot.png" target="_blank" >Plot of signal versus log2 fold change on insertion site flanking sequences</a><br>
<a href="$sessiondir/ta_densityplot.png " target="_blank" >Density Plot of log2 fold change with putative fold change cutoff(s)</a><br>
<a href="$sessiondir/ta_alloutputmerged.tsv" target="_blank" >Combined table of conditionally essential insertion site flanking sequences (.tsv)</a><br><br>
Add Control insertion site flanking sequences to database for comparison of different experiments (also required for MINOMICS visualization):
END
	$output=0;
	$output = add2db_link::add2db_html_link('motif','insertionsitecontrol_output_add2db.html',$sessiondir.'insertionsitecontrol_output_add2db.html',$workdir.'control_minomics.tsv','motif','Multiple','1','resource=notappl2','group1=none','group2=none');
	print OUTPUT " ".$output." \n";
	print OUTPUT <<END;
<br>Add Target insertion site flanking sequences to database for comparison of different experiments (also required for MINOMICS visualization):
END
	$output=0;
	$output = add2db_link::add2db_html_link('motif','insertionsitetarget_output_add2db.html',$sessiondir.'insertionsitetarget_output_add2db.html',$workdir.'target_minomics.tsv','motif','Multiple','1','resource=notappl3','group1=none','group2=none');
	print OUTPUT " ".$output." \n <hr>";
    } else {print OUTPUT "No conditionally essential insertion site analysis performed. Did you include target and control samples? <br><hr> \n";}
    
    if (-e "essentialgenes_alloutputmerged.tsv") {
	print OUTPUT <<END;
<b>Essential genes</b> <br>
<a href="$sessiondir/essentialgenes_PCA_non_normalized.png " target="_blank" >PCA plot of non normalized read counts on essential genes</a><br>
<a href="$sessiondir/essentialgenes_PCA_normalized.png" target="_blank" >PCA plot of normalized data read counts on essential genes</a><br>
<a href="$sessiondir/essentialgenes_plot.png " target="_blank" >Plot of signal versus log2 fold change on essential genes</a><br>
<a href="$sessiondir/essentialgenes_densityplot.png " target="_blank" >Density Plot of log2 fold change with putative fold change cutoff(s)</a><br>
<a href="$sessiondir/essentialgenes_alloutputmerged.tsv" target="_blank" >Combined table for essential genes (.tsv)</a><br><br>
Add essential gene data to database for comparison of different experiments (also required for MINOMICS visualization).<br>
END
	$output = 0;
	$output = add2db_link::add2db_html_link('madatabase','essentialgenes_output_add2db.html',$sessiondir.'essentialgenes_output_add2db.html',$workdir.'essentialgenes_minomics.tsv','madatabase','Multiple','1','resource=NA','group1=none','group2=none','basenumber=2');
	print OUTPUT " ".$output." \n";

	print OUTPUT <<END;
<hr>
<a href="http://bamics2.cmbi.ru.nl/websoftware/minomics" target="_blank" > Run MINOMICS on database data for visualization. Remember to add the data first. </a>
<br>
<br>
<hr>
END
    } else {print OUTPUT "Essential gene analysis failed. Check your input files \n";}

    if ($zip eq "Yes") {
	print OUTPUT <<END;
<br>
<hr>
<a href="$sessiondir/archive.zip" target="_blank" >Zipped file of output files excluding sequences and raw alignment output</a><br>    
END
    }
    close OUTPUT ;    
}#END generatehtml

#######################################################################################
# truncategenes
# This subroutine will chop of a small part of the 3' end of a gene to filter out matches with transposon insertions that have no effect on function
# it takes the log(1.1) of 10% of the length of the gene and removes that from the end. 
# This empirical value was shown to perform quite well with both long and short genes.

sub truncategenes {
    use POSIX;

    my $pttline;
    my $locusptt;
    my $startptt;
    my $stopptt;
    my $strandptt;
    my $pttvar4;
    my $pttvar5;
    my $pttvar6;
    my $pttvar7;
    my $pttlength;
    my $newstartptt;
    my $newstopptt;
    my @spliit;

    print "Applying gene truncation filter\n";
    open(PTT, 'genome.ptt') or general::error("Can't open file: genome.ptt");
    open (TRUNCATED, "> truncated.ptt");
    while (<PTT>){
        $pttline = $_ ;
        chomp $pttline;
	@spliit = split '\t', $pttline;
	$locusptt = $spliit[0] ;
	$startptt = $spliit[1] ;
	$stopptt  = $spliit[2] ;
	$strandptt = $spliit[3] ;
	$pttvar4 =  $spliit[4] ;
	$pttvar5 = $spliit[5];
	$pttvar6 = $spliit[6];
	$pttvar7 = $spliit[7];
	if ($locusptt eq "!Locus") {
	    print TRUNCATED "$locusptt\t$startptt\t$stopptt\t$strandptt\t$pttvar4\t$pttvar5\t$pttvar6\t$pttvar7\n";
	}
	next if ($locusptt eq "!Locus");
	$pttlength=$stopptt-$startptt;
	if ($pttlength <20) {
	    print TRUNCATED "$locusptt\t$startptt\t$stopptt\t$strandptt\t$pttvar4\t$pttvar5\t$pttvar6\t$pttvar7\n";
	}
	if ($pttlength >=20) {
	    if ($strandptt eq "-") {
		$newstartptt = $startptt;
		$newstartptt = floor($startptt + (log (0.05 * $pttlength))/(log(1.1)) );
		print TRUNCATED "$locusptt\t$newstartptt\t$stopptt\t$strandptt\t$pttvar4\t$pttvar5\t$pttvar6\t$pttvar7\n";
	    }
	    if ($strandptt eq "+") {
		$newstopptt = $stopptt;
		$newstopptt = floor($stopptt - (log (0.05 * $pttlength))/(log(1.1)) );
		print TRUNCATED "$locusptt\t$startptt\t$newstopptt\t$strandptt\t$pttvar4\t$pttvar5\t$pttvar6\t$pttvar7\n";
	    }
	}
    }
close PTT;
close TRUNCATED;
}#END truncategenes

#subroutine to count number of sequences in fasta file
sub countlines {
    my $filename = shift @_;
    my $lines = 0;
    my $buffer = 0;
    open(FILE, $filename) or die "Can't open `$filename': $!";
    while (<FILE>){ 
	if ($_ =~ m/\>/) {
	    $lines=$lines+1;}
	}
    close FILE;
    print "The number of sequences in $filename is $lines.\n";
}#END countlines

#subroutine to count number of sequences/lines in raw text file
sub countlines2 {

    my $filename = shift @_;
    my $lines = 0;
    my $buffer = 0;
    open(FILE, $filename) or die "Can't open `$filename': $!";
    while (sysread FILE, $buffer, 4096) {
	$lines += ($buffer =~ tr/\n//);
    }
    close FILE;
    my $sequences = $lines ;
    print "The number of sequences in $filename is $sequences.\n";
}#END countlines2


#countwiggle
sub countwiggle {
    my $wigglefile = shift @_;
    my $splet2;
    my $totallines = 0;
    my $totaldepth = 0;
    open(WIGGLEFILE, "< $wigglefile" ) or general::error("Can't open file: ".$wigglefile."");
    while (<WIGGLEFILE>){
	my $wigglefileline = $_ ;
	chomp $wigglefileline ;
	my @splet2 = split '\t', $wigglefileline;
	my $wigglefiledepth = $splet2[0];
	$totaldepth = $totaldepth + $wigglefiledepth;
	$totallines = $totallines + 1;
    }
    print "$totaldepth reads aligned to $totallines insertion site flanking sequences\n";
}#END countwiggle                                                                                                                                                