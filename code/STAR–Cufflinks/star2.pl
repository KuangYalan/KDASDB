#!/usr/bin/env perl
use strict;
use 5.010;
use warnings;

# use Data::Printer;
use File::Basename;
use IPC::Cmd qw[can_run run];
use File::Copy;
use Getopt::Long;

my $usage = <<USAGE;

SYSNOPSIS
star2.pl V1.1

This program is used to batch run star2.pl

star2.pl [options] genomeDir file|glob

 Options:
   -m|--2-pass        'none', 'multi' or 'per', default 'none'
   -e|--encode        ENCODE RNA-Seq Guide line
   -n|--mismatch      maximum mismatch number, default 6
   -s|--summarize     summarize only, default off
   -a|--args          argments used for star
   -g|--gtf           gtf file
   -o|--output        output folder for removed reads 
   -p|--progress      progress, default 16
   
Note: This program supported the STAR 2-pass mode

USAGE
my $gtf_file   = '';
my $encode     = 0;
my $help       = 0;
my $out_folder = dirname './';
my $mismatch   = 6;
my $summarize  = 0;
my $mode       = 'none';
my $progress   = 16;
my $args_file  = '';
die $usage
  unless GetOptions(
    "n|mismatch:s" => \$mismatch,
    "e|encode"     => \$encode,
    "a|args:s"     => \$args_file,
    "g|gtf:s"      => \$gtf_file,
    "s|summarize"  => \$summarize,
    "o|output:s"   => \$out_folder,
    "m|mode:s"     => \$mode,
    "p|progress=i" => \$progress,
    "h|help"       => \$help,
  );
##########################################################################################
#Checking the parameter infomation
##########################################################################################
die $usage if $help;
die "Not supported star 2-pass mode:$mode\n"
  if $mode ne 'per' && $mode ne "multi" && $mode ne "none";

#check the genome database
my $genomeDir = shift;
$genomeDir =~ s/\/$//;
die "Does not specifed STAR reference database\n" unless $genomeDir;
-s "$genomeDir/SA"
  || -s "$genomeDir/genomeParameters.txt"
  || die "Cannot detect STAR GenomeDir:$genomeDir\n";
my @fqs;

#check the genome fasta files
open my $genome_db_fh, "<", "$genomeDir/genomeParameters.txt";
my $genome_fa_file    = '';
my $genomeChrBinNbits = 18;
my $sjdbOverhang      = 100;
while (<$genome_db_fh>) {
    next if /^#/;
    chomp;
    my ( $key, $value ) = split /\t/;
    $genome_fa_file    = "$genomeDir/$value" if $key eq 'genomeFastaFiles';
    $genomeChrBinNbits = $value              if $key eq 'genomeChrBinNbits';
    $sjdbOverhang      = $value              if $key eq 'sjdbOverhang';
}

#check the fastq files;
foreach my $glob_str (@ARGV) {
    push @fqs, grep { -s $_ } glob $glob_str;
}
die $usage unless @fqs;

#check the running environment
can_run('STAR') or die 'STAR is not installed!';
can_run('zcat') or die 'zcat is not installed!';
if ($gtf_file) {
    die "$gtf_file does not exits\n" unless -s $gtf_file;
}

#mkdir
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;
my $log_folder = $out_folder . "/logs";
mkdir $log_folder unless -d $log_folder;

open my $results_fh, ">", "$out_folder/MappingStatistics.xls"
  or die "Canot create MappingStatistics.xls\n";

##########################################################################################
#Getting the mapping infomation
##########################################################################################

##########################################################################################
#Getting the paired-end or single-end data
##########################################################################################
my %fastqs = ();
my @suffx =
  qw (_1.fq.gz _2.fq.gz _1.fastq.gz _2.fastq.gz _1.fq _2.fq _1.fastq _2.fastq
  .1.fq.gz .2.fq.gz .1.fastq.gz .2.fastq.gz .1.fq .2.fq .1.fastq .2.fastq
  .fq.gz .fq .fastq.gz .fastq);
warn "Getting read pair information\n";
foreach my $fq ( sort @fqs ) {
    my $fullname = basename($fq);
    if ( $fq =~ /(.*)[\.|_]([1|2])\.(fq\.gz|fq|fastq\.gz|fastq)/ ) {
        my $basename = basename( $fq, @suffx );
        my $mate = $2;
        $fastqs{PE}{$basename}{$mate} = $fq;
        warn "$fullname\tPAIRED\t$mate\n";
    }
    elsif ( $fq =~ /(.*)\.(fq\.gz|fq|fastq\.gz|fastq)/ ) {
        my $basename = basename( $fq, @suffx );
        $fastqs{SE}{$basename} = $fq;
        warn "$fullname\tSINGLE\t1\n";

    }
    else {
        warn "Not supported format:$fq\n";
    }
}

my %reads_info = ();
foreach my $type ( sort keys %fastqs ) {

    if ( $type eq 'PE' ) {
        foreach my $basename ( sort keys %{ $fastqs{$type} } ) {
            my $mate1       = $fastqs{$type}{$basename}{1};
            my $mate2       = $fastqs{$type}{$basename}{2};
            my $fq_in_param = " -1 $mate1 -2 $mate2 ";

            $reads_info{$basename}{1} = $mate1;
            $reads_info{$basename}{2} = $mate2;

            $reads_info{$basename}{type}  = 'PE';
            $reads_info{$basename}{param} = $fq_in_param;

        }
    }
    else {
        foreach my $basename ( sort keys %{ $fastqs{$type} } ) {
            my $mate        = $fastqs{$type}{$basename};
            my $fq_in_param = " -U $mate ";
            $reads_info{$basename}{1}     = $mate;
            $reads_info{$basename}{type}  = 'SE';
            $reads_info{$basename}{param} = $fq_in_param;

        }
    }

}

#######################################################################
#construct general star papameter
#######################################################################
my $star_param = '';
$star_param .= "--sjdbGTFfile $gtf_file " if $gtf_file;
$star_param .= "--runThreadN $progress ";

$mismatch = 999 if $encode;
$star_param .= "--outFilterMismatchNmax $mismatch ";
$star_param .=
  "--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 "
  if $encode;
$star_param .=
"--alignSJDBoverhangMin 1 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 "
  if $encode;
$star_param .= "--alignIntronMax 1000000 --alignMatesGapMax 1000000 "
  if $encode;

# $star_param .= "--outReadsUnmapped Fastx ";
$star_param .=
  "--outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif";
if ($args_file) {
    open my $args_fh, "<", $args_file or die "cannot open $args_file\n";
    while (<$args_fh>) {
        chomp;
        $star_param .= "$_ ";
    }
}

my %fq_parameters = ();
foreach my $basename ( sort keys %reads_info ) {
    my $fq_in_param = "--readFilesIn";
    my $reads       = $reads_info{$basename}{type};
    if ( $reads eq "PE" ) {
        $fq_in_param .=
" $reads_info{$basename}{1}  $reads_info{$basename}{2} --outFilterIntronMotifs RemoveNoncanonical ";
    }
    else {
        $fq_in_param .= " $reads_info{$basename}{1} ";
    }
    $fq_in_param .= " --readFilesCommand zcat "
      if $reads_info{$basename}{1} =~ /\.gz$/;
    $fq_parameters{$basename} = $fq_in_param;
}

#######################################################################
#mkdir log files and output folders
#######################################################################
my $pass1_log_folder = "$log_folder/1-pass";
my $pass1_out_folder = "$out_folder/1-pass";
my $pass2_out_folder = "$out_folder/2-pass";
my $pass2_log_folder = "$log_folder/2-pass";

if ( $mode eq 'none' ) {
    warn "Staring star2 with 1-pass mode ...\n";
}
elsif ( $mode eq 'per' ) {
    warn "Staring star2 with per-sample 2-pass mode ...\n";
}
else {
    warn "Staring star2 with multi-sample 2-pass mode ...\n";
}
my $header = "Sample\tInput\tUnique\tMutliple\n";
print $results_fh $header;
my @pass1_junctions_db_files = ();
warn "########\n1-pass...\n########\n" if $mode eq 'multi';

foreach my $basename ( sort keys %fq_parameters ) {
    warn "$basename\n";

    #mkdir
    if ( $mode eq 'per' ) {
        mkdir $pass2_out_folder unless -d $pass2_out_folder;
        mkdir $pass2_log_folder unless -d $pass2_log_folder;
        my $sample_pass2_out_folder = "$pass2_out_folder/$basename";
        mkdir $sample_pass2_out_folder unless -d $sample_pass2_out_folder;

        my $pass1_command =
"STAR --genomeDir $genomeDir --twopassMode Basic $fq_parameters{$basename} $star_param --outFileNamePrefix $sample_pass2_out_folder/ >$pass2_log_folder/$basename.stdout.log 2> $pass2_log_folder/$basename.stderr.log";
        $pass1_command =~ s/\s+/ /g;
        warn "\t$pass1_command\n";
        run_command($pass1_command) unless $summarize;
        my %counts = summarize_star2("$sample_pass2_out_folder/Log.final.out");
        print $results_fh
          "$basename\t$counts{input}\t$counts{unique}\t$counts{muitiple}\n";
    }
    else {
        mkdir $pass1_out_folder unless -d $pass1_out_folder;
        mkdir $pass1_log_folder unless -d $pass1_log_folder;
        my $sample_pass1_out_folder = "$pass1_out_folder/$basename";
        mkdir $sample_pass1_out_folder unless -d $sample_pass1_out_folder;

        my $pass1_command =
"STAR --genomeDir $genomeDir $fq_parameters{$basename} $star_param --outFileNamePrefix $sample_pass1_out_folder/ >$pass1_log_folder/$basename.stdout.log 2> $pass1_log_folder/$basename.stderr.log";
        $pass1_command =~ s/\s+/ /g;
        warn "\t$pass1_command\n";
        run_command($pass1_command) unless $summarize;
        if ( $mode eq 'none' ) {
            my %counts =
              summarize_star2("$sample_pass1_out_folder/Log.final.out");
            print $results_fh
              "$basename\t$counts{input}\t$counts{unique}\t$counts{muitiple}\n";
        }
        push @pass1_junctions_db_files, "$sample_pass1_out_folder/SJ.out.tab";
    }
}

##################################################################
################multi-sample 2-pass align
##################################################################

if ( $mode eq 'multi' ) {
    warn "########\n2-pass...\n########\n";
    mkdir $pass2_out_folder unless -d $pass2_out_folder;
    mkdir $pass2_log_folder unless -d $pass2_log_folder;

    #get the juncion db files
    my $pass1_junc_db_file_str = join " ", @pass1_junctions_db_files;

    #get the juncion db files
    foreach my $basename ( sort keys %fq_parameters ) {

        #get the input fastq parameter
        warn "$basename\n";

        my $sample_pass2_out_folder = "$pass2_out_folder/$basename";
        mkdir $sample_pass2_out_folder unless -d $sample_pass2_out_folder;

        my $pass2_command =
"STAR --genomeDir $genomeDir $fq_parameters{$basename} --sjdbFileChrStartEnd $pass1_junc_db_file_str $star_param --outFileNamePrefix $sample_pass2_out_folder/ >$pass2_log_folder/$basename.stdout.log 2> $pass2_log_folder/$basename.stderr.log";
        $pass2_command =~ s/\s+/ /g;
        warn "\t$pass2_command\n";
        run_command($pass2_command) unless $summarize;
        my %counts = summarize_star2("$sample_pass2_out_folder/Log.final.out");
        print $results_fh
          "$basename\t$counts{input}\t$counts{unique}\t$counts{muitiple}\n";
    }
}

##################################################################
################reformat the files
##################################################################
my $bam_folder = "$out_folder/Bams";
mkdir $bam_folder unless -e $bam_folder;
my $junction_folder = "$out_folder/Junctions";
mkdir $junction_folder unless -e $junction_folder;

foreach my $basename ( sort keys %fq_parameters ) {
    if ( $mode eq 'none' ) {
        move( "$pass1_out_folder/$basename/Aligned.sortedByCoord.out.bam",
            "$bam_folder/$basename.bam" );
        move(
            "$pass1_out_folder/$basename/SJ.out.tab",
            "$junction_folder/$basename.tab"
        );

        # move ($pass1_out_folder, $details_folder)
    }
    elsif ( $mode eq 'per' ) {
        move( "$pass2_out_folder/$basename/Aligned.sortedByCoord.out.bam",
            "$bam_folder/$basename.bam" );
        move(
            "$pass2_out_folder/$basename/SJ.out.tab",
            "$junction_folder/$basename.tab"
        );
    }
    elsif ( $mode eq 'multi' ) {
        unlink "$pass1_out_folder/$basename/Aligned.sortedByCoord.out.bam";
        move( "$pass2_out_folder/$basename/Aligned.sortedByCoord.out.bam",
            "$bam_folder/$basename.bam" );
        move(
            "$pass2_out_folder/$basename/SJ.out.tab",
            "$junction_folder/$basename.tab"
        );
    }

}

sub make_pass_dir {
    my $lable = shift;

}

sub run_command {
    my $command = shift;
    my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
      run( command => $command, verbose => 0 );
    if ($success) {

        # warn "\tDone!\n";
    }
    else {
        warn "Something went wrong:\n@$stderr_buf";
    }
}

sub summarize_star2 {
    my $star2_file = shift;
    open my $fh, "<", $star2_file or die "cannot open $star2_file\n$!";
    my %data;
    while (<$fh>) {
        chomp;
        next unless $_;
        s/^\s+//;
        s/\s+\|//;
        my ( $key, $value ) = split /\t/;
        $data{$key} = $value;
    }
    close $fh;
    my %counts = (
        input    => $data{'Number of input reads'},
        unique   => $data{'Uniquely mapped reads number'},
        muitiple => $data{'Number of reads mapped to multiple loci'},
    );
    return %counts;
}

