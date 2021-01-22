#!/usr/local/perl

#0: filename

$gff_file=$ARGV[0]."-FGS.gff";
$faa_file=$ARGV[0]."-FGS.faa";
$ffn_file=$ARGV[0]."-FGS.ffn";
$fna_file=$ARGV[0].".fna";


print "$gff_file";


$check_gff="Done";
%contig=();

open (SEQ,"$fna_file");
@tmp0=split(/\t/,$SEQ);
$col2=$tmp0[1];
if($col2==''){
$col2='FGS';
}
while (<SEQ>){
    chomp $_;
    if ($_=~/\>(.+)/){
        $h=$1;      # Modified on Jan 10 2018
        @varh = split  ' ', $h, 2;
        $head=$varh[0];
        $contig{$head}="";
    }
    else{
        $contig{$head}=$contig{$head}.$_;
    }

}
close(SEQ);


%cds=();
%ncbiname=();
%ncbigene=();
open (FILE,"$gff_file");
while (<FILE>){
    $aline = $_;
    if ($aline=~/^\#/){
    }else{
	#@tmp=split(/\s+/,$aline);
	#YY March 2019
	@tmp=split(/\t/,$aline);
	if ($tmp[2]=~/CDS/){
	    if (length($tmp[0]>40)){
		$tmp[0]=substr $tmp[0],0,40;
	    }
	    if ($tmp[3]=~/\d+/ && $tmp[4]=~/\d+/ && $tmp[6]=~/[+|-]/){
		$value=$tmp[3]."_".$tmp[4]."_".$tmp[6];
	    }else{
		$check_gff="error:gff format error";
		print "$check_gff\n";
		exit;
	    }
	    @tmp2=split(/\t/, $aline);
 	    if ($tmp2[8] =~ /Name=([^\;]+)/) {
		$ncbiname0 = $1;
	    } else {
		$ncbiname0 = "unk";
	    }
 	    if ($tmp2[8] =~ /gene=([^\;]+)/) {
		$ncbigene0 = $1;
	    } else {
		$ncbigene0 = "unk";
	    }
	    if (exists $cds{$tmp[0]}){
		$cds{$tmp[0]}=$cds{$tmp[0]}.";".$value;
		$ncbiname{$tmp[0]}=$ncbiname{$tmp[0]}.";".$ncbiname0;
		$ncbigene{$tmp[0]}=$ncbigene{$tmp[0]}.";".$ncbigene0;
	    }else{
		$cds{$tmp[0]}=$value;
		$ncbiname{$tmp[0]}=$ncbiname0;
		$ncbigene{$tmp[0]}=$ncbigene0;
	    }
	    
	}
    }
}
close(FILE);


$gff_size=keys %cds;
$fasta_size=keys %contig;
if ($gff_size==0){
    $check_gff="Warning:gff file did not contain CDS information, no gene?";
    print "$check_gff\n";
  #  exit;
}

if ($fasta_size>$gff_size){
    $check_gff="Warning:gff file did not contain all input sequence information.";
    print "$check_gff\n";
#    exit;
}

$new_gff=$gff_file.".new";
open (FAA,">$faa_file");
open (NEWGFF,">$new_gff");
open(FFN,">$ffn_file");
for $contig_id (keys %contig){
    $contig_flag=0;
    for $cds_id(keys %cds){
	$check=index($contig_id,$cds_id);
	if ($check!=-1){
	    $contig_flag=1;
	    @each_cds=split(/\;/,$cds{$cds_id});
	    @each_names=split(/\;/, $ncbiname{$cds_id});
	    @each_genes=split(/\;/, $ncbigene{$cds_id});
	    for ($i=0;$i<=@each_cds-1;$i++){
		@tmp=split(/\_/,$each_cds[$i]);
		$contig{$contig_id}=~ s/\r//g;          # Added on Jan 10 2018
                $gene=substr $contig{$contig_id}.$contig{$contig_id},($tmp[0]-1),($tmp[1]-$tmp[0]+1);
		if ($tmp[2]=~/\-/){
		    $gene=&comp($gene);
		}
		$protein = "";
		for ($j = 0; $j < (length($gene) - 2); $j = $j + 3){
		    $codon=substr($gene, $j, 3);
                    $aa = &codon2aa($codon);
                    #if ($j==0){            # Added on Jan 10 2018 
                     #   $aa='M';
                    #}
		    $aa = &codon2aa($codon);
		    if (($aa eq "*") && ($j < (length($gene) - 3))){
			$aa="U";
		    }
		    $protein = $protein.$aa;
		}
		$protein=~s/\*//g;
		if (length($protein)==0){
		    $check_gff="error: gene position error in provided gff file; skipped";
		    print "$check_gff\n";
		    next; 
		}
		$protein_id=$contig_id."_".$tmp[0]."_".$tmp[1]."_".$tmp[2];
		$old_id = $each_names[$i];
		$old_gene = $each_genes[$i];
		print NEWGFF "$contig_id\tFGS\tCDS\t$tmp[0]\t$tmp[1]\t.\t$tmp[2]\t.\tID=$protein_id\;Name=$protein_id\;Ori=$old_id\;gene=$old_gene\;Parent\n";
		print FFN ">$protein_id\n";
                print FAA ">$protein_id\n";
		my @gene_70 = unpack("(A70)*", $gene);
                for ($line=0;$line<=@gene_70-1;$line++){
                    print FFN "$gene_70[$line]\n";
                }
                my @protein_70 = unpack("(A70)*", $protein);
		#for ($line=0;$line<=@protein_70-1;$line++){
		    #print FAA "$protein_70[$line]\n";
		#}
		print FAA "$protein\n";
	    }
	    last;
	}
    }
    if ($contig_flag==0){
	$check_gff="Warning:no gene information for $contig_id";
	print "$check_gff\n";
#	exit;
    }
}
close(FNA);
close(FAA);

#use File::Copy qw(copy);
rename $new_gff, $gff_file;

print "$check_gff\n";


sub codon2aa{
    my($codon)=@_;
    $codon=uc $codon;
my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
    if(exists $g{$codon})
    {
	return $g{$codon};
    }
    else
    {
	return "X";
    }
}

sub comp{
    my ($str_nu)=@_;

    @nuclo=split(//,$str_nu);
    for ($n=0;$n<=@nuclo-1;$n++){
	if ($nuclo[$n]=~/T/){
            $nuclo[$n]=~ s/T/A/;
	}
        elsif ($nuclo[$n]=~/A/){
            $nuclo[$n]=~ s/A/T/;
        }
        elsif ($nuclo[$n]=~/C/){
            $nuclo[$n]=~ s/C/G/;
        }
        elsif ($nuclo[$n]=~/G/){
            $nuclo[$n]=~ s/G/C/;
        }
    }
    @re=reverse(@nuclo);
    $re_str=join("",@re);
    return ($re_str);
}
close(NEWGFF);
