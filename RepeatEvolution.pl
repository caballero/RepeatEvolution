#!/usr/bin/perl

=head1 NAME

RepeatEvolution.pl

=head1 DESCRIPTION

Synthetic evolution of repeat consensi.

=head1 USAGE

RepeatEvolution.pl [OPTIONS]

    Parameter       Description                     Value       Default
    -f --fasta      Fasta file                      File
    -o --out        Output file                     File
    -m --matrix     Scoring matrix                  File
    -c --config     Configuration file              File
    -s --saveseqs   Save intermediate sequences
    
    -v --verbose    Verbose mode
    -h --help       Print this screen

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
my $matrix    = undef;

# Main variables
my $our_version = 0.1;
my %seq;
my %conf;
my %matrix;
my %ins_size_p;
my %del_size_p;

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'f|fasta=s'        => \$fasta,
    'c|config=s'       => \$config,
    'o|out=s'          => \$out,
    'm|matrix=s'       => \$matrix,
    's|saveseqs'       => \$save_seqs
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);
printVersion() if (defined $version);

pod2usage(-verbose => 2) if !(defined $fasta);
pod2usage(-verbose => 2) if !(defined $config);
pod2usage(-verbose => 2) if !(defined $out);
pod2usage(-verbose => 2) if !(defined $matrix);

loadConfig($config);
loadFasta($fasta);
loadMatrix($matrix);

my $mut_rate   = $conf{    'mutation_rate'};
my $max_expand = $conf{    'max_expansion'};
my $dead_lim   = $conf{   'dead_threshold'};
my $last_gen   = $conf{'total_generations'};
my $report     = $conf{'report_generation'};

my $snv_lim    = $conf{         'snv_freq'};
my $del_lim    = $conf{         'del_freq'} + $snv_lim;
my $ins_lim    = $conf{         'ins_freq'} + $del_lim;

# MAIN LOOP
for (my $gen = 1; $gen <= $last_gen; $gen++) {
    warn "GENERATION $gen\n" if (defined $verbose and $gen % $report == 0);
    foreach my $id (keys %seq) {
        next if ($seq{$id}{'isDead'} == 1); # sequence is dead
        # do we want to replicate?
        if ($seq{$id}{'exp'} >= rand) {
            # how many copies?
            my $num_copies = int(rand $max_expand);
            $num_copies = 1 if ($num_copies < 1);
            warn "REPLICATING $id in $num_copies\n" if (defined $verbose);
            for (my $j = 1; $j <= $num_copies; $j++) {
                my $new_id = "$id-$j";
                $seq{$new_id}{'div'}       =      $seq{$id}{'div'};
                $seq{$new_id}{'exp'}       =      $seq{$id}{'exp'};
                $seq{$new_id}{'mut'}       =      $seq{$id}{'mut'};
                $seq{$new_id}{'isDead'}    =   $seq{$id}{'isDead'};
                @{ $seq{$new_id}{'seq'} }  = @{ $seq{$id}{'seq'} };
            }
        }
        mutate($id, $mut_rate);
        $seq{$id}{'isDead'} = 1 if ($seq{$id}{'div'} > $dead_lim);
    }
    saveSeqs("$fasta.G$gen.fa") if (defined $save_seqs and $gen % $report == 0);
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
            $id  =   $1;
            $id .= '-0';
            $seq{$id}{'div'}      =   0;
            $seq{$id}{'mut'}      =   0;
            $seq{$id}{'exp'}      = $conf{'prob_expansion'};
            $seq{$id}{'isDead'}   =   0;
            @{ $seq{$id}{'seq'} } =  ();
        }
        else {
            my $seq = uc $_;
               $seq = validateSeq($seq) if ($seq =~ m/[^ACGT]/);
            my @seq = split(//, $seq);
            push @{ $seq{$id}{'seq'} }, @seq;
        }
    }
    close F;
    warn "found $cnt sequences\n" if (defined $verbose);
}

sub loadConfig {
    my ($file) = @_;
    warn "reading configuration from $file\n" if (defined $verbose);
    open F, "$file" or die "cannot open $file\n";
    while (<F>) {
        chomp;
        next if (m/^#/ or m/^\n/); # skip comments and empty lines
        my ($var, $val) = split (/\s*:=\s*/, $_);
        $conf{$var} = $val;
    }
    close F;
    
    my @ins_p      = split (/,/, $conf{'ins_dist'});
    foreach my $ins (@ins_p) {
        my ($s, $f) = split (/=/, $ins);
        $ins_size_p{$s} = $f;
    }
    
    my @del_p      = split (/,/, $conf{'del_dist'});
    foreach my $del (@del_p) {
        my ($s, $f) = split (/=/, $del);
        $del_size_p{$s} = $f;
    }
    
}

sub loadMatrix {
    my ($file) = @_;
    warn "reading score matrix from $file\n" if (defined $verbose);
    open F, "$file" or die "cannot open $file\n";
    my $pos = 0;
    while (<F>) {
        chomp;
        next if (m/^#/);
        my ($a, $c, $g, $t, $i, $d) = split (/\t/, $_);
        $matrix{$pos}{'A'} = $a;
        $matrix{$pos}{'C'} = $c;
        $matrix{$pos}{'G'} = $g;
        $matrix{$pos}{'T'} = $t;
        $matrix{$pos}{'I'} = $i;
        $matrix{$pos}{'D'} = $d;
        $matrix{$pos}{'N'} =  0; # neutral mutation
        $pos++;
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
        my $div = sprintf("%.2f", 100 * $seq{$id}{'div'});
        my $seq = formatSeq(join "", @{ $seq{$id}{'seq'} });
        print F ">$id | $div\%\n$seq";
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
    my $myrate  = rand($rate);
    my @myseq   = @{ $seq{$id}{'seq'} };
    my $num_mut = int($myrate * (1 + $#myseq));
    
    for (my $i = 1; $i <= $num_mut; $i++) {
        my $pos = int(rand @myseq);
        my $dice = rand;
        if    ($dice < $snv_lim) { 
            my $b = addSNV($id, $pos);
            $seq{$id}{'exp'} += $matrix{$pos}{$b};
        }
        elsif ($dice < $del_lim) { 
            addDel($id, $pos);
            $seq{$id}{'exp'} += $matrix{$pos}{'D'};
        }
        elsif ($dice < $ins_lim) { 
            addIns($id, $pos);
            $seq{$id}{'exp'} += $matrix{$pos}{'I'};
        }   
    }
    
    my $mut   = $seq{$id}{'mut'};
    my $div   = sprintf("%.4f", $mut / ($#myseq + 1));
    $seq{$id}{'div'} = $div;
}

sub addSNV {
    my ($id, $pos) = @_;
    my $b = $seq{$id}{'seq'}[$pos];
    if (length $b == 1) {
        # regular SNV
        my @n = qw/ a c g t/;
        if    ($b =~ m/a/i) { @n = qw/c g t/; }
        elsif ($b =~ m/c/i) { @n = qw/a g t/; }
        elsif ($b =~ m/g/i) { @n = qw/a c t/; }
        elsif ($b =~ m/t/i) { @n = qw/a c g/; }
        my $n = $n[ int(rand @n) ];
        $seq{$id}{'seq'}[$pos] = $n;
        $seq{$id}{'mut'}++;
        return uc $n;
    }
    elsif (length $b > 1) {
        # mutation in an insertion point
        my $pos2 = int(rand(length $b));
        my $a    = substr($b, $pos2, 1);
        my @n = qw/ a c g t/;
        if    ($a =~ m/a/i) { @n = qw/c g t/; }
        elsif ($a =~ m/c/i) { @n = qw/a g t/; }
        elsif ($a =~ m/g/i) { @n = qw/a c t/; }
        elsif ($a =~ m/t/i) { @n = qw/a c g/; }
        substr($b, $pos2, 1) = $n[ int(rand @n) ];
        $seq{$id}{'seq'}[$pos] = $b;
        return 'N';
    }
    else {
        # do nothing is a deletion point
        return 'N';
    }
}

sub addDel {
    my ($id, $pos) = @_;
    my $max_size = 1;
    my $dice = rand;
    my $sum  = 0;
    foreach my $del (keys %del_size_p) {
        $sum     += $del_size_p{$del};
        $max_size = $del if ($sum < $dice);
    }
    my $size = int(rand $max_size);
    $size = 1 if ($size < 1);
    for (my $i = 1; $i <= $size; $i++) {
        $seq{$id}{'seq'}[$pos] = '';
    }
    $seq{$id}{'mut'} += $size;
}

sub addIns {
    my ($id, $pos) = @_;
    my $max_size = 1;
    my $dice = rand;
    my $sum  = 0;
    foreach my $ins (keys %ins_size_p) {
        $sum     += $ins_size_p{$ins};
        $max_size = $ins if ($sum < $dice);
    }
    my $size = int(rand $max_size);
    $size = 1 if ($size < 1);
    $seq{$id}{'seq'}[$pos] .= newSeq($size);
    $seq{$id}{'mut'} += $size;
}

sub newSeq {
    my ($len) = @_;
    my $seq   = '';
    my @n     = qw/a a t t g c/;
    for (my $i = 1; $i <= $len; $i++) {
        $seq .= $n[ int(rand @n) ];
    }
    return $seq;
}
