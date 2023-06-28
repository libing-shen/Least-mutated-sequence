#!/usr/bin/perl -w

# Program: MutationSiteFinder v.1.0
# Find mutation site in DNA alignment and output a matrix
# Libing Shen
# International Human Phenome Institutes (Shanghai)
# Email contact <shenlibing@ihup.org.cn>
# June 2023

use strict;
use warnings;
use Data::Dumper;

my $input_file = shift @ARGV;

my %pos_homolog;

my $pos_ref =\%pos_homolog;

($input_file, $pos_ref) = read_fasta_and_create_hash($input_file, $pos_ref);

my @mutation_site;

my $mut_ref =\@mutation_site;

($pos_ref, $mut_ref) = find_mutation_site($pos_ref, $mut_ref);

output_matrix_file($mut_ref, $pos_ref);

exit;

sub read_fasta_and_create_hash {
    
    my ($input_file, $pos_ref) = @_;
    
    $/ = ">";
    
    my $INPUT_FILE_HANDLE;
    
    open ($INPUT_FILE_HANDLE, $input_file) || die $!;
    
    my $first_item = <$INPUT_FILE_HANDLE>;
    
    while (my $item = <$INPUT_FILE_HANDLE>){
        
        chomp $item;
        
        my @lines = split /\n/, $item;
        
        my $header = shift @lines;
        
        $header =~ s/\r//g;
        
        my $sequence;
        
        foreach my $seq_bit (@lines){
            
            $sequence .= $seq_bit;
        }
        
        $sequence =~ s/\r//g;
        
        my @sequence = split /|/, $sequence; 
        
        for (my $i = 0; $i < scalar(@sequence); $i++){
            
            my %array_elment = (
                "header" => "$header",  
                "nucleotide" => "$sequence[$i]",   
            );
            
            my $hash_ref =\%array_elment;
           
            push (@{${$pos_ref}{$i}}, $hash_ref);
        }
    }

    close $INPUT_FILE_HANDLE;
    
    return ($input_file, $pos_ref);
}

sub find_mutation_site {
    
    my ($pos_ref, $mut_ref) = @_;
    
    foreach my $key (sort {$a <=> $b} keys %{$pos_ref}){
        
        my %mutation;
        
        foreach my $element (@{$pos_homolog{$key}}){
            
            my $nucleotide = $element->{nucleotide};
            
            $mutation{$nucleotide} = 0;
        }
        
        my @nucleotide = keys (%mutation);
        
        if (scalar(@nucleotide) != 1){
            
            push (@{$mut_ref}, $key);
        }
    }
    
    return ($pos_ref, $mut_ref);
}

sub output_matrix_file {
    
    my ($mut_ref, $pos_ref) = @_;
    
    my @mut_site = @{$mut_ref};
    my %pos_ref = %{$pos_ref};
    
    my $first = $mut_site[0];
    
    my $output_file = "mutation_site_matrix.txt";
    
    my $OUTPUT_FILE_HANDLE;
    
    open ($OUTPUT_FILE_HANDLE, ">$output_file") || die $!;
    
    foreach my $element (@mut_site){
        
        if ($element == $first){
            
            print $OUTPUT_FILE_HANDLE "position";
            
            foreach my $item (@{$pos_ref{$element}}){
                
                print $OUTPUT_FILE_HANDLE "\t".$item->{header};
            }
            
            print $OUTPUT_FILE_HANDLE "\n";
            
            my $pos = $element+1;
            
            print $OUTPUT_FILE_HANDLE $pos;
            
            foreach my $item (@{$pos_ref{$element}}){
                
                print $OUTPUT_FILE_HANDLE "\t".$item->{nucleotide};
            }
            
            print $OUTPUT_FILE_HANDLE "\n";
            
        } else {
            
            my $pos = $element+1;
            
            print $OUTPUT_FILE_HANDLE $pos;
            
            foreach my $item (@{$pos_ref{$element}}){
                
                print $OUTPUT_FILE_HANDLE "\t".$item->{nucleotide};
            }
            
            print $OUTPUT_FILE_HANDLE "\n";   
        }
    }
    
    close $OUTPUT_FILE_HANDLE;
}
