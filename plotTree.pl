#!/usr/bin/perl

=head1 NAME

plotTree.pl

=head1 DESCRIPTION

=head1 USAGE

=head1 EXAMPLES

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
my $help     = undef;         # Print help
my $verbose  = undef;         # Verbose mode
my $version  = undef;         # Version call flag

# Main variables
my $our_version = 0.1;        # Script version number

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose
) or pod2usage(-verbose => 2);
    
pod2usage(-verbose => 2) if (defined $help);
printVersion() if (defined $version);

my %nodes;
my %pairs;
my %nchildren;
my %names;
my $max_child = 0;
my %levels;

while (<>) {
    chomp;
    next unless (m/^>/);
    s/>//;
    my ($name, $div) = split (/ \| /, $_);
    unless ($name =~ m/-0$/) {
        my $parent = $name;
        $parent =~ s/-\d+$//;
        $pairs{$parent}{$name} = 1;
        $nchildren{$parent}++;
        $max_child = $nchildren{$parent} if ($max_child < $nchildren{$parent});
    }
    $div =~ s/\%//;
    $nodes{$name} = $div;
    my $fname = $name;
    $fname =~ s/-//g;
    $names{$name} = $fname;
}


print "graph G {\n    node [shape=circle, style=filled, fixedsize=true, colorscheme=rdylbu10];\n";

foreach my $node (sort keys %nodes) {
    my $chd = 0;
    $chd = $nchildren{$node} if (defined $nchildren{$node});
    next if ($chd == 0);
    
    my @div = ();
    push @div, $nodes{$node};
    my $min_div = $nodes{$node};
    my $max_div = $nodes{$node};
    my $sum_div = $nodes{$node};
    my $num_div = 1;
        foreach my $child (keys %{ $pairs{$node} }) {
        push @div, $nodes{$child};
        $sum_div += $nodes{$child};
        $min_div  = $nodes{$child} if ($min_div > $nodes{$child});
        $max_div  = $nodes{$child} if ($max_div < $nodes{$child});
        $num_div ++;
    }
    my $mean_div = sprintf("%.2f", $sum_div / $num_div);
    #my $divRange = 0.6;
    #my $color = sprintf("%0.2f", (($mean_div - $min_div) * 0.6) / $divRange);
    my $color = int(3 * $mean_div / 10) + 1;
    my $name = $names{$node};
    $chd++;
    my $size = sprintf("%.2f", 2 * sqrt($chd / $max_child));
    print "    $name [label=\"$node\\ndiv=$mean_div\%\\nN=$chd\", height=$size, width=$size, color=$color];\n";
    #$levels{int($mean_div / 2)} .= "    $name;\n";
}

#foreach my $level (keys %levels) {
#    print "    subgraph div$level {\n        rank=same;\n        $levels{$level} }\n";
#}

foreach my $parent (sort keys %pairs) {
    my $pname = $names{$parent};
    foreach my $child (sort keys %{ $pairs{$parent} }) {
        next unless (defined $nchildren{$child});
        my $cname = $names{$child};
        print "    $pname -- $cname;\n";
    }
}
print "}";

###################################
####   S U B R O U T I N E S   ####
###################################

# printVersion => return version number
sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}

