#!/usr/local/perl

#0: fasta 
#1: output fasta
#2: summary 

#Running this script, will remove the header, and only keeps the accession number
#For instance NZ_FIZP01000001.1 Campylobacter geochelonis strain RC20, whole genome shotgun sequence 
#will become NZ_FIZP01000001.1

print "$ARGV[0]";

$seq_check_out="Done";
open (CONT,"$ARGV[0]");
%contig_seq=();
%contig_des=();
$flag=0;
$n=1;
$short = "Query_1";
$full = "Query_1";
while (<CONT>){
    chomp $_;
    if ($_=~/^\>(.+)/){
	$flag++;
	@head=split(/\s+/,$1);
	$short = shift @head;
 	$full = join(" ", @head);
	if (length($short)==0){
	    $seq_check_out="error1: format error";
	    print "$seq_check_out\n";
	    exit;
	}
	if (length($short)>40){
	    $short=substr $short,0,40;
	}
	#$short=$short."_".$n;
	$short=$short; #to check YY June 12, 2018
	$n++;
        $contig_seq{$short}="";
	$contig_des{$short}=$full;
    }
    else{
	$_=uc($_);
	if ($_!~/^[A|T|G|C|N]+$/){
	    @letter=split(//,$_);
	    for ($i=0;$i<=@letter-1;$i++){
		 if ($letter[$i]=~/[A|T|G|C|N]/){
		     $contig_seq{$short}=$contig_seq{$short}.$letter[$i]; 
		 }elsif ($letter[$i]=~/[M|R|W|V|H|D|S|Y|B|K|X]/){
		     $contig_seq{$short}=$contig_seq{$short}."N";  
		 }else{
		     $contig_seq{$short}=$contig_seq{$short}.$letter[$i];
		 }
	    }
	}else{
	    $contig_seq{$short}=$contig_seq{$short}.$_;
	}
    }
}
close (CONT);
#add strictfasta, some users submit sequences that have no fasta header >
$strictfasta = 0;
if ($strictfasta==1 && $flag==0){
    $seq_check_out="error2: not fasta format";
    print "$seq_check_out\n";
    exit;
}else{
    $tot = 0;
    for $key (keys %contig_seq){
	if (length($contig_seq{$key})==0){
	    $seq_check_out="error3: format error";
	    print "$seq_check_out\n";
	    exit;
	}
	$contig_seq{$key}=~s/[^\w]//g;
        if ($contig_seq{$key}!~/^[A|T|G|C|N|X]+$/){
	    $seq_check_out="error4: You need to input nucleotide sequence";
            print "$seq_check_out\n";
            exit;
        }
	$tot=$tot+1;
    }
    if($tot==0) {
    	$seq_check_out="error2: not fasta format";
    	print "$seq_check_out\n";
    	exit;
    }
}

print "$seq_check_out\n";

open (NEW,">$ARGV[1]");
open (SUM,">$ARGV[2]");
for $key (keys %contig_seq){
   print NEW ">$key\n";
   my @seq = unpack("(A70)*", $contig_seq{$key});
   for ($i=0;$i<=@seq-1;$i++){
       print NEW "$seq[$i]\n";
   }
   $len=length($contig_seq{$key});
   $note=$contig_des{$key};
   print SUM "$key\t$len\t$note\n";
}
close(NEW);
close(SUM);
