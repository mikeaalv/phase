#!/usr/bin/perl
open(INFILE, "/Users/mikeaalv/Documents/miRNA/Phase_RNA_results/wp/root/T_AM_sly_cds.sam");
open(OUTFILE, ">/Users/mikeaalv/Documents/miRNA/Phase_RNA_results/wp/root/T_AM_sly_cds_pl")|| die "Cannot open the newfile: $!\n";;
while (<INFILE>) {
        @a = split(" ");
        my $l=length($a[5]);
        my $long=($a[4])+$l;
        my $i = $a[4] +1;
        print OUTFILE "$a[0]\t$a[3]\t$i\t$long\t$a[2]\t$a[1]\n";
}
exit;
