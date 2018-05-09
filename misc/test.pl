#!/usr/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use DBI;


=head1 NAME

=head1 SYNOPSIS

Options:

	-help       brief help message
	-man        full documentation

=head1 DESCRIPTION

=cut

# parse command line options
my $help;
my $man;
my $result = GetOptions (	"help"	=> \$help,
							"man"	=> \$man);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);



my @ary = DBI->available_drivers();
print join("\n", @ary), "\n";

