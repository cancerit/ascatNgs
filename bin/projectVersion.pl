#!/usr/bin/env perl

=head
Innanely simple script to get the tagged version
Can be ported to any project if put in the same directory structure
=cut

use strict;
use Cwd qw(abs_path);
use File::Basename;
my $path = dirname(abs_path($0)).'/../version.properties';
my $line = `grep currentVersion $path`;
$line =~ s/^currentVersion=//;
print $line;
