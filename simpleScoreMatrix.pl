#!/usr/bin/perl

=head1 NAME

simpleScoreMatrix.pl

=head1 DESCRIPTION

Create a basic score matrix to use in RepeatEvolution.pl.

=head1 USAGE

simpleScoreMatrix.pl [OPTIONS]

    Parameter       Description             Value       Default
    -f --fasta      Fasta file              File        STDIN
    -o --out        Output file             File        STDOUT
    -m --match      Score for match         [0,1]       0.01
    -x --mismatch   Score for mistmatch     [0,1]       -0.01
    -i --insertion  Score for insertion     [0,1]       -0.01
    -d --deletion   Score for deletion      [0,1]       -0.01

=head1 EXAMPLES

    1. Simple usage:
    perl simpleScoreMatrix.pl < FASTA > SCORE.matrix
    
    2. Calling files
    perl simpleScoreMatrix.pl -f FASTA -o SCORE.matrix
    
    3. Changing scores
    perl simpleScoreMatrix.pl -f FASTA -o SCORE.matrix -m 0.1 -x 0 -i 0 -d 0
   

=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2012

=head1 CONTACT

jcaballero@systemsbiology.org

=head1 LICENSE

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with code.  If not, see <http://www.gnu.org/licenses/>.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Default parameters
my $help      = undef;         # Print help
my $verbose   = undef;         # Verbose mode
my $version   = undef;         # Version call flag
my $fasta     = undef;
my $out       = undef;
my $match     =  0.01;
my $mismatch  = -0.01;
my $insertion = -0.01;
my $deletion  = -0.01;

# Main variables
my $our_version = 0.1;        # Script version number
my @dna = qw/A C G T/;;


# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'version'          => \$version,
    'f|fasta:s'        => \$fasta,
    'o|out:s'          => \$out,
    'm|match:f'        => \$match,
    'x|mismatch:f'     => \$mismatch,
    'i|insertion:f'    => \$insertion,
    'd|deletion:f'     => \$deletion 
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);
printVersion() if (defined $version);

if (defined $fasta) {
    open STDIN, "$fasta" or die "cannot open file $fasta\n";
}

if (defined $out) {
    open STDOUT, ">$out" or die "cannot open file $out\n";
}

while (<>) {
    chomp;
    if (m/>/) {
        s/>//;
        print "#SEQ: $_\n";
        print "#"; print join "\t", @dna, 'I', 'D'; print "\n";
    }
    else {
        my @seq = split (//, uc $_);
        foreach my $b (@seq) {
            my @val = ();
            if ($b =~ m/N/) {
                @val = (0, 0, 0, 0, 0, 0);
            }
            else {
                $b = decode($b) if ($b =~ m/[^ACGT]/);
                foreach my $dna (@dna) {
                    if ($b =~ m/$dna/) { 
                        push @val, $match; 
                    }
                    else { 
                        push @val, $mismatch; 
                    }
                }
                push @val, $insertion; 
                push @val, $deletion;
            }
            print join "\t", @val;
            print "\n";
        }
    }   
}

###################################
####   S U B R O U T I N E S   ####
###################################

# printVersion => return version number
sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}

sub decode {
    my ($n) = @_;
    if    ($n eq 'S') { $n =  'CG'; }
    elsif ($n eq 'R') { $n =  'AG'; }
    elsif ($n eq 'Y') { $n =  'CT'; }
    elsif ($n eq 'W') { $n =  'AT'; }
    elsif ($n eq 'K') { $n =  'GT'; }
    elsif ($n eq 'M') { $n =  'AC'; }
    elsif ($n eq 'B') { $n = 'CGT'; }
    elsif ($n eq 'D') { $n = 'AGT'; }
    elsif ($n eq 'H') { $n = 'ACT'; }
    elsif ($n eq 'V') { $n = 'ACG'; }
    return $n;
}    
