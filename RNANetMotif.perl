#!/usr/bin/perl

use strict;
use warnings;

# command line arguments
# 1. input file directory
# 2. k
my $inputdir = $ARGV[0];
my $k = $ARGV[1];

my @inputdirArray = split(/\//, $inputdir);
my $fafile = pop(@inputdirArray);
my $workingdir = join("\/", @inputdirArray) . "\/";

my @nameArray=split(/\./,$fafile);
my $RBP=$nameArray[0];


my $E_file=$workingdir."E_profile_".$RBP."_u1.txt";
my $H_file=$workingdir."H_profile_".$RBP."_u1.txt";
my $I_file=$workingdir."I_profile_".$RBP."_u1.txt";
my $M_file=$workingdir."M_profile_".$RBP."_u1.txt";
my $output=$workingdir."unpaired_probability_".$RBP."_u1.txt";
   
system("E_RNAplfold -W 100 -L 100 -u 1 <$inputdir >$E_file");
system("H_RNAplfold -W 100 -L 100 -u 1 <$inputdir >$H_file");
system("I_RNAplfold -W 100 -L 100 -u 1 <$inputdir >$I_file");
system("M_RNAplfold -W 100 -L 100 -u 1 <$inputdir >$M_file");
system("python combine_letter_profiles.py $E_file $H_file $I_file $M_file 1 $output");
system("rm ".$workingdir."*profile*txt");

system("cd $workingdir && RNAplfold -W 100 -L 100 -c 0.5 -u 1 <$inputdir && cd \~-");
system("rm ".$workingdir."*_lunp");

system("rename '(' '_' ".$workingdir."*");
system("rename ')' '' ".$workingdir."*");

my @psfiles=glob("$workingdir*.ps");

foreach my $psfilePath(@psfiles)
{
   my @temp=split(/\//,$psfilePath);
   my $psfile =pop(@temp);
   my @nameArray=split(/_/,$psfile);
   my $outfile;
   if( $#nameArray>1)
   {
   	$outfile = join('_', @nameArray[0..$#nameArray-2])."_edges";
   }
   else
   {
        $outfile = join('_', @nameArray[0..$#nameArray-1])."_edges";
   }
   #my $outfile = join('_', @nameArray[0..$#nameArray-2])."_edges";
   my $outfilePath = $workingdir.$outfile;
   system("cd $workingdir && "."grep -E '[0-9] ubox' $psfile >$outfile");
}

system("rm ".$workingdir."*ps");

system("python graphk_kmer.py ".$workingdir." ".$RBP." ".$k);

system("python convert_faseq_edges_to_dataframe_new.py"." ".$workingdir." ".$inputdir);
system("python convert_faseq_basepair_to_dataframe_new.py"." ".$workingdir." ".$inputdir);

system("rm ".$workingdir."*edges");
system("python generate_nodeattri_transtographfeature.py ".$k." ".$workingdir." ".$fafile);


system("python vdm3distance_findcliques_motifcsv.py ".$workingdir." ".$RBP." ".$k." 0 1 6 1 ".$fafile);

