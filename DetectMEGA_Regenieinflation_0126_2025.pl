#Created it on 0410_2024
#Collect all INDs that maybe involved into Regenie GWAS analysis
#Exclude the problem IDs due to age and flagged bad quality
#AllINDs_ForRegenieGWAS Source
#1)Individuals from the MEGA dataset included in the MEGA imputed data.
#2)QCed Individuals from the New MEGA EUR imputed data
#3)Exclude individuals with their age >100
#3)Remove the INDs with bad quality which was checked by new MEGA
#GWAS regenie GWAS result fille
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
#1 727233 rs151190501 A G 0.978579 0.659311 85589 ADD 0.0897587 0.0933994 0.923559 0.472961 NA

#!/usr/local/perl

use List::Util qw( min max );
use POSIX qw/ceil/;
use Text::CSV; 
use List::MoreUtils qw(first_index);

use strict;
use warnings;

sub ReadAllChiSQ_Pvalues_inRegenieGWASResultFile
{
    my ($RegenieBasedGWASFile,$chr)=@_;
    my %rsid_Chisq_pvalue=();
    my $LineNumber=0;
    open (INPUT, "gunzip -c $RegenieBasedGWASFile|") or die "gunzip $RegenieBasedGWASFile: $!";;
     while(<INPUT>)
     {
         my @cols = split(/\s++/,$_,2);
	 if($LineNumber==0)
	 {	     	 
            $LineNumber++;
            next;
	 }
       if($cols[0] == $chr)	 
       {
        @cols = split(/\s++/,$_);
        $rsid_Chisq_pvalue{$cols[2]}="$cols[11]:$cols[12]";
        $LineNumber++;
       }
       last if ($cols[0] >$chr);
     }
     close INPUT;
     my $size=keys %rsid_Chisq_pvalue;
     return (\%rsid_Chisq_pvalue);
}


sub FilterDifferentPvalue_cHISQ_inRegenieGWASFile
{
    my ($RegenieBasedGWASFile,$rsid_ChiSQ_pvalueRef,$chr)=@_;
    my %rsid_ChiSQ_pvalue_OneFile=%{$rsid_ChiSQ_pvalueRef};
    my %rsid_ChiSQ_pvalue_Filtered=();
    my %rsid_A1Freq_INFO_Filtered=();
    my $LineNumber=0;
    open (INPUT, "gunzip -c $RegenieBasedGWASFile|") or die "gunzip $RegenieBasedGWASFile: $!";;
     while(<INPUT>)
     {
	  my @cols = split(/\s++/,$_,2);    
         if($LineNumber==0)
         {
            $LineNumber++;
            next;
         }
	if($cols[0] == $chr)
       {
         @cols = split(/\s++/,$_);
	 my $ChiSQ_pvalue="$cols[11]:$cols[12]"; 
	 if($rsid_ChiSQ_pvalue_OneFile{$cols[2]} ne $ChiSQ_pvalue)
	 {
         $rsid_ChiSQ_pvalue_Filtered{$cols[2]}="$cols[11]:$cols[12]";
         $rsid_A1Freq_INFO_Filtered{$cols[2]}="$cols[5]\t$cols[6]";
        }
       }
       last if ($cols[0] >$chr);
       $LineNumber++;
       #last if($LineNumber>1000);
     }
     close INPUT;
     my $RawSize=keys %rsid_ChiSQ_pvalue_OneFile;
     my $Filteredsize=keys %rsid_ChiSQ_pvalue_Filtered;
     print "$chr: GWAS kept variants:size: $RawSize vs $Filteredsize\n";
     return (\%rsid_ChiSQ_pvalue_Filtered,\%rsid_A1Freq_INFO_Filtered);
}


#GWAS regenie GWAS result fille
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
#1 727233 rs151190501 A G 0.978579 0.659311 85589 ADD 0.0897587 0.0933994 0.923559 0.472961 NA
## Make the Manhattan plot on the gwasResults dataset
#manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P" )

#Create it on 01_24_2025
#FilterVariantByMAF and convert it to Mahatton-plot 
sub FilterRegenieGWASResultFileByMAF
{
    my ($RegenieBasedGWASFile,$MAFThreshold,$FilteredGWASFile)=@_;
    my $LineNumber=0;
     my $UpperAF=1.0-$MAFThreshold;
     open OUT, ">$FilteredGWASFile" or die "Can't open Output file:$FilteredGWASFile!";
     print OUT "CHR\tBP\tSNP\tP\n";
     open (INPUT, "gunzip -c $RegenieBasedGWASFile|") or die "gunzip $RegenieBasedGWASFile: $!";;
     while(<INPUT>)
     {
         my @cols = split(/\s++/,$_);
         if($LineNumber==0)
         {
            $LineNumber++;
	    print OUT "$cols[0]\t$cols[1]\t$cols[2]\tP\n";
            next;
         }
	 #next if($cols[12] eq "NA");
	 
         if($MAFThreshold<$cols[5] && $cols[5]<$UpperAF)
	 {
          my $pvalue = 10 ** ((-1)* $cols[12]);
            if($pvalue>0.05)	
	    {
	       my $formatted_pvalue = sprintf("%.4f", $pvalue);
               print OUT "$cols[0]\t$cols[1]\t$cols[2]\t$formatted_pvalue\n";	       
	    }  
	    else
	    {
	      print OUT "$cols[0]\t$cols[1]\t$cols[2]\t$pvalue\n";
            }
	  $LineNumber++;
	 }
	 #last if($LineNumber>100);
     }
     close INPUT;
    close OUT;
}
#FIRTH/EUR_batch_16_EUR_153.3.regenie.gz   FIRTH/EUR_batch_3_EUR_041.4.regenie.gz  SPA/MEGA_v61_spaEUR_568.regenie.gz
#FIRTH/EUR_batch_170_EUR_568.regenie.gz    SPA/MEGA_v61_spaEUR_041.4.regenie.gz    SPA/MEGA_v61_spaEUR_687.1.regenie.gz
#FIRTH/EUR_batch_196_EUR_687.1.regenie.gz  SPA/MEGA_v61_spaEUR_153.3.regenie.gz
#FIRTH/EUR_batch_20_EUR_172.22.regenie.gz  SPA/MEGA_v61_spaEUR_172.22.regenie.gz

 my @PhecodeList=("041.4","153.3","172.22","568","687.1");
my @FirthBasedList=(3,16,20,170,196);

=cut;
#Construct two different methods based GWAS p value joint file(FIRTH and SPA)

for my $iPhecode (0..$#PhecodeList) {
 my $FirthBasedGWASFile= "FIRTH/EUR_batch_$FirthBasedList[$iPhecode]_EUR_$PhecodeList[$iPhecode].regenie.gz";
 my $SPABasedGWASFile= "SPA/MEGA_v61_spaEUR_$PhecodeList[$iPhecode].regenie.gz";
#Construct two different methods based GWAS p value different 
  my $InflatedVariantBasedInfoList="EUR_InflatedBasedFilterListInfo_$PhecodeList[$iPhecode].txt";
  open OUT, ">$InflatedVariantBasedInfoList" or die "Can't open Output file:$InflatedVariantBasedInfoList!";
 print OUT "rsid\tChiSQ_Firth\tLog10P_Firth\tChiSQ_SPA\t";
 print OUT "lOG10P_SPA\tA1Freq\tINFO\n";

  my @ChrList=(1 .. 23);
 for my $chr (1..$#ChrList) {
  my ($rsid_ChiSQ_pvalue_FirthRef)=ReadAllChiSQ_Pvalues_inRegenieGWASResultFile($FirthBasedGWASFile,$chr);
  my ($rsid_ChiSQ_pvalue_FilteredRef,$rsid_A1Freq_INFO_FilteredRef)=FilterDifferentPvalue_cHISQ_inRegenieGWASFile($SPABasedGWASFile,$rsid_ChiSQ_pvalue_FirthRef,$chr);
  my %rsid_ChiSQ_pvalue_Firth=%{$rsid_ChiSQ_pvalue_FirthRef};
  my %rsid_ChiSQ_pvalue_Filtered=%{$rsid_ChiSQ_pvalue_FilteredRef};
  my %rsid_A1Freq_INFO_Filtered=%{$rsid_A1Freq_INFO_FilteredRef};

  foreach my $rsid (sort keys %rsid_ChiSQ_pvalue_Filtered)
 {
  my @FirthBasedChiSQ_Pvalue=split(/:/,$rsid_ChiSQ_pvalue_Firth{$rsid});
  my @ChiSQ_pvalue_Filtered=split(/:/,$rsid_ChiSQ_pvalue_Filtered{$rsid});
  print OUT "$rsid\t$FirthBasedChiSQ_Pvalue[0]\t$FirthBasedChiSQ_Pvalue[1]\t$ChiSQ_pvalue_Filtered[0]\t";
  print OUT "$ChiSQ_pvalue_Filtered[1]\t$rsid_A1Freq_INFO_Filtered{$rsid}\n";
  }
 }
 close OUT;
 }
=cut;

#MAF_Filtering :use MAF to filter GWAS summary statistics
my @MAFThresholds=(0.01);
for my $iMAF (0..$#MAFThresholds) {
	for my $iPhecode (0..$#PhecodeList) {
 		my $FirthBasedGWASFile= "FIRTH/EUR_batch_$FirthBasedList[$iPhecode]_EUR_$PhecodeList[$iPhecode].regenie.gz";
 	        my $FirthFilteredGWASFile="FIRTH/MAF_Filtering/EUR_$PhecodeList[$iPhecode].regenie_MAFThreshold$MAFThresholds[$iMAF].txt";
 		FilterRegenieGWASResultFileByMAF($FirthBasedGWASFile,$MAFThresholds[$iMAF],$FirthFilteredGWASFile);
	}
}
