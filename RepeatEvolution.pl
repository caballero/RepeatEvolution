#!/usr/bin/perl

=head1 NAME

RepeatEvolution.pl

=head1 DESCRIPTION

Synthetic evolution of repeat consensi.

=head1 USAGE

RepeatEvolution.pl -f FASTA -c CONFIG -o OUPUT [-s]

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
my $help      = undef;
my $verbose   = undef;
my $version   = undef;
my $fasta     = undef;
my $config    = undef;
my $out       = undef;
my $save_seqs = undef;

# Main variables
my $our_version = 0.1;
my %seq;
my %conf;

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'f|fasta:s'        => \$fasta,
    'c|config:s'       => \$config,
    'o|out:s'          => \$out,
    's|saveseqs'       => \$save_seqs
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);
printVersion() if (defined $version);

pod2usage(-verbose => 2) if !(defined $fasta);
pod2usage(-verbose => 2) if !(defined $config);
pod2usage(-verbose => 2) if !(defined $out);

loadFasta($fasta);
loadConfig($config);

my $mut_rate   = $conf{    'mutation_rate'};
my $p_expand   = $conf{   'prob_expansion'};
my $max_expand = $conf{    'max_expansion'};
my $dead_lim   = $conf{   'dead_threshold'};
my $last_gen   = $conf{'total_generations'};
my $report     = $conf{'report_generation'};

# MAIN LOOP
for (my $gen = 1; $gen <= $last_gen; $gen++) {
    warn "GENERATION $gen\n" if (defined $verbose and ($gen % $report == 0) );
    foreach my $id (keys %seq) {
        next if ($seq{$id}{'div'} < $dead_lim); # sequence is dead;
        # do we want to replicate?
        if ($p_expand >= rand) {
            # how many copies?
            my $num_copies = int(rand $max_expand);
            warn "REPLICATING $id in $num_copies\n" if (defined $verbose);
            for (my $j = 1; $j <= $num_copies; $j++) {
                my $new_id = "$id-$j";
                $seq{$new_id}{'div'}  = $seq{$id}{'div'};
                $seq{$new_id}{'seq'}  = $seq{$id}{'seq'};
            }
        }
        
        mutate($id, $mut_rate);
        
        saveSeqs("$fasta.G$gen.fa") if (defined $save_seqs);
    }
}

# final sequences after simulation
saveSeqs($out);

###################################
####   S U B R O U T I N E S   ####
###################################

sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}

sub defineFH {
    my ($f) = @_;
    my $fh  = $f;
    $fh = "gunzip  -c $f | " if ($f =~ m/gz$/);
    $fh = "bunzip2 -c $f | " if ($f =~ m/bz2$/);
    return $fh;
}

sub loadFasta {
    my ($file) = @_;
    warn "reading sequences from $file\n" if (defined $verbose);
    my $fh = defineFH($file);
    open F, "$fh" or die "cannot open $file\n";
    my $id;
    my $cnt = 0;
    while (<F>) {
        chomp;
        if (m/>(.+)/) {
            $cnt++;
            $id = $1;
            $seq{$id}{'div'}  = 0.0;
            $seq{$id}{'seq'}  =  '';
        }
        else {
            my $seq = validateSeq(uc $_);
            $seq{$id}{'seq'} .= $seq;
        }
    }
    close F;
    warn "found $cnt sequences\n" if (defined $verbose);
}

sub loadConfig {
    my ($file) = @_;
    open F, "$file" or die "cannot open $file\n";
    while (<F>) {
        chomp;
        next if (m/^[#\n]/); # skip comments and empty lines
        my ($var, $val) = split (/\s*:=\s*/, $_);
        $conf{$var} = $val;
    }
    close F;
}

sub validateSeq {
    my ($s) = @_;
    return $s unless ($s =~ m/[^ACGT]/);
    $s =~ s/U/T/;
    my @s = split (//, $s);
    for (my $i = 0; $i <= length $#s; $i++) {
        next if ($s[$i] =~ m/[ACGT]/);
        my @b = qw/A C G T/;
        if    ($s[$i] eq 'S') { @b = qw/C G/; }
        elsif ($s[$i] eq 'R') { @b = qw/A G/; }
        elsif ($s[$i] eq 'Y') { @b = qw/C T/; }
        elsif ($s[$i] eq 'W') { @b = qw/A T/; }
        elsif ($s[$i] eq 'K') { @b = qw/G T/; }
        elsif ($s[$i] eq 'M') { @b = qw/A C/; }
        elsif ($s[$i] eq 'B') { @b = qw/C G T/; }
        elsif ($s[$i] eq 'D') { @b = qw/A G T/; }
        elsif ($s[$i] eq 'H') { @b = qw/A C T/; }
        elsif ($s[$i] eq 'V') { @b = qw/A C G/; }
        
        $s[$i] = $b[int(rand @b)];
    }
    return join "", @s;
}

sub saveSeqs {
    my $file = shift @_;
    warn "saving sequence in $file\n" if (defined $verbose);
    open F, ">$file" or die "cannot open $file\n";
    foreach my $id (keys %seq) {
        my $div = $seq{$id}{'div'};
        my $seq = formatSeq($seq{$id}{'seq'});
        print F ">$id | $div\n$seq";
    }
    close F;
}

sub formatSeq {
    my ($seq) = @_;
    my $for   = '';
    while ($seq) {
        $for .= substr($seq, 0, 80);
        $for .= "\n";
        substr($seq, 0, 80) = '';
    }
    return $for;
}

sub mutate {
    my ($id, $rate) = @_;
    
    
}
