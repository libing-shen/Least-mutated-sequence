#!/usr/bin/perl -w

# Program: MutationSiteCounter v.1.0
# Count point mutation sites in matrix and output a sorted result
# Libing Shen
# International Human Phenome Institutes (Shanghai)
# Email contact <shenlibing@ihup.org.cn>
# June 2023

use strict;
use warnings;
use Data::Dumper;

my $input_file  = shift @ARGV; 

#my $input_file = "mutation_site_matrix.txt";

my @sample;
my @position;
my %sequence;

my $sample_ref =\@sample;
my $position_ref =\@position;
my $sequence_ref =\%sequence;

($input_file, $sample_ref, $position_ref, $sequence_ref) = read_matrix_file_and_process ($input_file, $sample_ref, $position_ref, $sequence_ref);

my %reduced_seq1;
my %reduced_seq2;

my $reduced_seq1 =\%reduced_seq1;
my $reduced_seq2 =\%reduced_seq2;

($sample_ref, $sequence_ref, $reduced_seq1, $reduced_seq2) = merge_identical_sequence ($sample_ref, $sequence_ref, $reduced_seq1, $reduced_seq2);

my @parsimony;

my $parsimony =\@parsimony;

($parsimony, $reduced_seq1, $reduced_seq2) = calculate_least_mutated_seq ($parsimony, $reduced_seq1, $reduced_seq2);

my $output_file = "mutation_count_result.txt";;

my $OUTPUT_FILE_HANDLE;

open ($OUTPUT_FILE_HANDLE, ">$output_file") || die $!;

print $OUTPUT_FILE_HANDLE "point_mutations\ttransitions\ttransversions\tTiTv_ratio\tcases\n";

foreach my $element (@parsimony){
    
    my $titv_ratio = 0;
    
    if ($element->{tv} != 0){
        
        $titv_ratio = $element->{ti}/$element->{tv};
    }
    
    $titv_ratio =  sprintf("%.2f", $titv_ratio);
    
    print $OUTPUT_FILE_HANDLE $element->{point}."\t".$element->{ti}."\t".$element->{tv}."\t".$titv_ratio."\t".$element->{name}."\n";
}

exit;

sub read_matrix_file_and_process {
    
    my ($input_file, $sample_ref, $position_ref, $sequence_ref) = @_;
    
    my $INPUT_FILE_HANDLE;
    
    open ($INPUT_FILE_HANDLE, $input_file) || die $!;
    
    print STDERR "Start processing \n";

    my $counter = 0;
    
    while (my $line = <$INPUT_FILE_HANDLE>){
        
        chomp $line;
        
        if ($line =~ /^position/){
            
            my @words = split /\t/, $line;
            
            shift @words;
            
            @{$sample_ref} = @words;
            
        } else{
            
            my @words = split /\t/, $line;
            
            my $pos = shift @words;
            
            push (@{$position_ref}, $pos);
            
            for (my $i = 0; $i < scalar(@words); $i++){
                
                if (!exists $sequence{$i}){
                    
                    push (@{${$sequence_ref}{$i}}, $words[$i])
                    
                } else{
                    
                    push (@{${$sequence_ref}{$i}}, $words[$i])
                }
            }
            
            $counter++;
            
            print STDERR "process line ".$counter."\n";
        }
    }
    
    close $INPUT_FILE_HANDLE;
    
    return ($input_file, $sample_ref, $position_ref, $sequence_ref);
}


sub merge_identical_sequence {
    
    my ($sample_ref, $sequence_ref, $reduced_seq1, $reduced_seq2) = @_;
    
    print STDERR "Start merging \n";
    
    my @merge;

    for (my $i = 0; $i < scalar(@{$sample_ref}); $i++){
        
        my $array_ref =\@{${$sequence_ref}{$i}};
        
        my $string = join ("", @{${$sequence_ref}{$i}});
        
        my %array_elment =(
            "name" => $sample[$i],
            "seq" => $array_ref,
            "string" => $string,
        );
        
        my $hash_ref =\%array_elment;
        
        push (@merge, $hash_ref);
    }
    
    print STDERR "Start reducing \n";

    my %reduced;

    foreach my $element (@merge){
        
        if (!exists $reduced{$element->{string}}){
            
            push (@{$reduced{$element->{string}}}, $element->{name});
            
        } else{
            
            push (@{$reduced{$element->{string}}}, $element->{name});
        }
    }
    
    foreach my $key (keys %reduced){
        
        my $name = join (",", @{$reduced{$key}});
        
        ${$reduced_seq1}{$name} = $key;
    }
    
    %{$reduced_seq2} = %{$reduced_seq1};
    
    return ($sample_ref, $sequence_ref, $reduced_seq1, $reduced_seq2);
}

sub calculate_least_mutated_seq {
    
    my ($parsimony, $reduced_seq1, $reduced_seq2) = @_;
    
    print STDERR "Start calculating \n";
    
    my $cal_counter = 0;
    
    foreach my $key1 (keys %{$reduced_seq1}){
        
        my $point_mutation = 0;
        my $transition = 0;
        my $transversion = 0;
        
        my @seq2 = split /|/, ${$reduced_seq1}{$key1};
        
        foreach my $key2 (keys %{$reduced_seq2}){
            
            if ($key1 ne $key2){
                
                my @seq3 = split /|/, ${$reduced_seq2}{$key2};
                
                for (my $i = 0; $i < scalar(@seq2); $i++){
                    
                    if ($seq2[$i] ne $seq3[$i]){
                        
                        my @array;
                        
                        push (@array, $seq2[$i]);
                        push (@array, $seq3[$i]);
                        
                        @array = sort{$a cmp $b} @array;
                        
                        my $match = join("", @array);
                        
                        if ($match eq "AG" || $match eq "CT"){
                            
                            $transition++;
                            
                        } elsif ($match eq "AT" || $match eq "AC" || $match eq "CG" || $match eq "GT"){
                            
                            $transversion++;
                            
                        } else{
                            
                            print Dumper ($match);
                        }
                        
                        $point_mutation++;
                    }
                }
            }
        }
        
        my %array_elment =(
            "name" => $key1,
            "point" => $point_mutation,
            "ti" => $transition,
            "tv" => $transversion,
        );
        
        my $hash_ref =\%array_elment;
        
        push (@{$parsimony}, $hash_ref);
        
        $cal_counter++;
        
        print STDERR "Finishing sequence ".$cal_counter." ".localtime."\n";
    }
    
    @{$parsimony} = sort {$a->{point} <=> $b->{point}} @{$parsimony};

    return ($parsimony, $reduced_seq1, $reduced_seq2);
}
