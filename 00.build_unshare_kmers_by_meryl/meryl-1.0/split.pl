#!/bin/perl

my @inputFiles;
while (scalar(@ARGV)) {
    my $arg = shift @ARGV;
    push @inputFiles, $arg;
}

my $fileNumber    = "001";
my $fileLength    = 0;
my $fileLengthMax = 10000000000;
open(OUT, "| gzip -1c > reads-$fileNumber.fasta.gz");
print STDERR "--   --> 'reads-$fileNumber.fasta.gz'.\n";

foreach my $file (@inputFiles) {

    print STDERR "--   <-- '$file'.\n";
    open(INP, "< $file");
    open(INP, "gzip  -dc $file |")   if ($file =~ m/\.gz$/);
    open(INP, "bzip2 -dc $file |")   if ($file =~ m/\.ba2$/);
    open(INP, "xz    -dc $file |")   if ($file =~ m/\.xz$/);

    while (!eof(INP)) {
        $_ = <INP>;

        #  If looks like FASTA or FASTQ name, make new sequence
        if ((m/^\@/) ||
            (m/^>/)) {
            print OUT ">\n";
            next;
        }

        #  If looks like FASTQ quality values, skip them
        if (m/^\+/) {
            $_ = <INP>;
            next;
        }

        #  Otherwise, assume it's sequence to keep, print it
        print OUT $_;
        $fileLength += length($_) - 1;
        if ($fileLength > $fileLengthMax) {
            close(OUT);
            $fileNumber++;
            $fileLength = 0;
            print STDERR "--   --> 'reads-$fileNumber.fasta.gz'.\n";
            open(OUT, "| gzip -1c > reads-$fileNumber.fasta.gz");
        }
    }
    close(INP);
}
close(OUT);
