use strict;
use warnings;
use Getopt::Long;

my $suffix = ".maskShort.fas";
my $filelist = "fastalist";
my $minlen = 9; # sequences shorter than this will be masked
my $mindist = 6; # if flanked on both sides with more than this number of gaps

GetOptions ("filelist=s" => \$filelist,   
            "suffix=s"   => \$suffix,    
            "minlen=i"   => \$minlen,    
            "mindist=i"   => \$mindist)    
or die("Error in command line arguments\n");

my @files;
open my $fh, '<', $filelist or die "cant open list of files to check: $!\n";
while(<$fh>){
	chomp;
	push @files, $_;
}
close($fh);

print "<maskshortseqs.pl> Read " , scalar @files , " file names to check from file $filelist\n";
print "Settings: minlen=", $minlen, ", dist=", $mindist,", suffix=",$suffix,"\n";

my $count = 0;
my $bpmasked = 0;
foreach my $file (@files){
  my $outfile = $file.$suffix;
  if (!-e $file){
    print "Warning:cant open for reading file $file, skipping\n";
    next;
  }
  open my $fhi, '<', $file or die "Warning:cant open for reading file $file: $!\n";
  
  $count++;
  open my $fho ,'>', $outfile or die "cant open for writing file $outfile:$!\n"; 
  while(<$fhi>){
    if (/^>/){
      print $fho $_;
    }
    else{
    	$_ =~ s/ +//g;
      my $nleftgaps = 0;
      my $nrightgaps = 0;
      my $ninside = 0;
      my $startinpos = -1;
      my @temp = split(//,$_);
      for my $i (0..$#temp){
        if($temp[$i] eq "-" || $temp[$i] eq "n" || $temp[$i] eq "N"){
          if($ninside == 0){
            $nleftgaps++;
          }
          else{
            $nrightgaps++;
          }
        }
        else{
          if($nleftgaps != 0){
              if($startinpos == -1){
                $startinpos = $i;
              }
              $ninside++;
            if($nrightgaps != 0){
              if($ninside < $minlen && $nleftgaps > $mindist && $nrightgaps > $mindist){
                for my $j ($startinpos..$startinpos+$ninside-2){
                  $temp[$j] = 'n';
                  $bpmasked++;
                }
              }
              $nleftgaps = $nrightgaps;
              $startinpos = $i;
              $nrightgaps = 0;
              $ninside = 1;
            }
          }
        }
      }
      my $final = join(@temp,"");
      print $fho @temp;

    }
  }
  close($fhi);
  close($fho);
} 

print "Number of files checked: $count, number of bp masked: $bpmasked\n";