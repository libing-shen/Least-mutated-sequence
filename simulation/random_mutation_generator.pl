#!/usr/bin/perl -w
use strict;
use Data::Dumper;

$/ = ">";

my $input_ancestral_sequence = "ancestral_sequence.fa";

my $INPUT_FILE_HANDLE;

open ($INPUT_FILE_HANDLE, "$input_ancestral_sequence") || die $!;

my $first_item = <$INPUT_FILE_HANDLE>;

my @mutated_sequence;

my $ancestor_sequnece;

while (my $item = <$INPUT_FILE_HANDLE>){
    
    chomp $item;
    
    my @line = split /\n/, $item;
    
    my $head = shift @line;
    
    my $sequence;
    
    foreach my $line (@line){
        
        $sequence .= $line;
    }
    
    $ancestor_sequnece = substr($sequence, 0, 100);
    
    my @nucleotides = ("A","C","G","T");
    
    srand(time|$$);
    
    for (my $j = 0; $j< 100; $j++){
        
	my $origin_seq = $ancestor_sequnece;
        
        my $iteration = $j + 1;
        
        for (my $i = 0; $i< 200; $i++){
            
            my $position = int(rand(length($origin_seq)));
            
            my $mut;
            
            do {
                
                $mut = $nucleotides[rand(@nucleotides)];
                
            } until ($mut ne substr($origin_seq, $position, 1));
            
            substr($origin_seq, $position, 1, $mut);
            
            my $mut_num = $i+1;
            
            my $name = "iter:".$iteration."_mut:".$mut_num;
            
            my %array_element = (
                "name" => "$name",
                "sequence" => "$origin_seq",
            );
            
            my $hash_ref =\%array_element;
            
            push (@mutated_sequence, $hash_ref);
        }
    }
}

close $INPUT_FILE_HANDLE;

my @output_array;

for (my $i = 0; $i < 10 ;$i++){
    
    my $rand_element = $mutated_sequence[int(rand(@mutated_sequence))];
    
    push (@output_array, $rand_element);
}

my $sim_num = "0001";

my $output_file = "simulation_".$sim_num.".fa";

my $OUTPUT_FILE_HANDLE;

open ($OUTPUT_FILE_HANDLE, ">$output_file") || die $!;

foreach my $element (@output_array){
    
    my $name = $element->{name};
    
    print $OUTPUT_FILE_HANDLE ">AY274119.3_".$name."\n";
    
    my $sequence = $element->{sequence};
    
    for ( my $i = 0; $i < length($sequence ); $i += 70) {
        
        print $OUTPUT_FILE_HANDLE substr($sequence, $i, 70), "\n";
    }
}

close $OUTPUT_FILE_HANDLE;

exit;
