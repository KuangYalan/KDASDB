#!/usr/bin/env perl
use strict;
use 5.010;
use warnings;

# use Data::Printer;
use File::Basename;
use Cwd 'abs_path';
use IPC::Cmd qw[can_run run];
use Getopt::Long;

my $usage = <<USAGE;
SYSNOPSIS
stringtie.pl [options] file|glob

 Options:
   -s|--bam-flag           default NA, can be setted as'accepted_hits'
   -g|--genome-guided      genome reference gtf file for RABT, optional
   -m|--mask-file          genome sequence fasta file, optional
   -f|--min-isfm-fraction  multiple aligned read count correct, default 0.1
   -c|--min-coverage       minimum read coverage, default: 2.5 
   -a|--args               argments used for stringtie
   -o|--output             output folder
   -p|--progress           progress, default 16
USAGE
my $out_folder    = dirname './';
my $progress      = 16;
my $genome_guided = '';
my $bam_flag      = '';

my $min_isoform_fraction = 0.1;
my $min_coverage         = 2.5;
my $args_file            = '';
die $usage
  unless GetOptions(
    "a|args:s"                 => \$args_file,
    "g|genome-guided:s"        => \$genome_guided,
    "s|bam-flag:s"             => \$bam_flag,
    "f|min-isoform-fraction:f" => \$min_isoform_fraction,
    "c|min-coverage"           => \$min_coverage,
    "o|output:s"               => \$out_folder,
    "p|progress=i"             => \$progress,
  );

file_check($genome_guided) if $genome_guided;

can_run('stringtie') or die 'stringtie is not installed!';
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;
my $log_folder = $out_folder . "/logs";
mkdir $log_folder unless -d $log_folder;

my @bams = ();
foreach my $glob_str (@ARGV) {

    #pick up directory
    my @bam_tmp = ();
    if ($bam_flag) {
        @bam_tmp = map { abs_path("$_/$bam_flag.bam") }
          grep { -d $_ and -s "$_/$bam_flag.bam" } glob $glob_str;
    }
    else {
        @bam_tmp = map { abs_path($_) }
          grep { -s $_ } glob $glob_str;
    }

    #pick up bams
    push @bams, @bam_tmp if @bam_tmp;
}

die "No bams were found, please check the glob string\n"
  unless @bams;

my $cufflinks_param = "-p $progress -f $min_isoform_fraction -c $min_coverage ";
$cufflinks_param .= "-G $genome_guided " if $genome_guided;

if ($args_file) {
    open my $args_fh, "<", $args_file or die "cannot open $args_file\n";
    while (<$args_fh>) {
        chomp;
        $cufflinks_param .= "$_ ";
    }
}

warn "Staring cufflinks...\n";
foreach my $bam ( sort @bams ) {
    my @paths = split /\//, $bam;

    my $basename = $bam;
    if ( $bam_flag && $bam =~ /$bam_flag.bam/ ) {
        $basename = $paths[-2];
    }
    else {
        $basename = basename( $bam, ".uniq.bam", "_uniq.bam", ".bam" );
        warn "$bam does not match bam flag, sample name set as:$basename \n";
    }

    my $stdout = "$basename.stdout.log";
    my $stderr = "$basename.stderr.log";
    my $command =
"stringtie -o $out_folder/$basename.gtf -l $basename $cufflinks_param $bam >$log_folder/$stdout 2> $log_folder/$stderr";

    warn "\t$command\n";
    my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
      run( command => $command, verbose => 0 );
    if ($success) {
        warn "\tDone!\n";
    }
    else {
        my @stderrs = @$stderr_buf;
        warn "Something went wrong:\n@stderrs";
    }

}

sub file_check {
    my $file = shift;
    die "File $file does not exists\n"
      unless -e $file;
    die "File:$file is empty\n"
      unless -s $file;
}
