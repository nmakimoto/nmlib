#! /usr/bin/perl
# Wrapper script of gtest_main

use Getopt::Long;
use strict;
use warnings;

my %OPT;
GetOptions(
  help => sub{print "Usage: $0 [--filter=PATTERN] [--list] [--repeat=COUNT] [--shuffle]\n"; exit 0},
  "filter=s" => \$OPT{filter},
  "list!"    => \$OPT{list},
  "repeat=i" => \$OPT{repeat},
  "shuffle!" => \$OPT{shuffle},
)  or  die;

my @cmd=("./gtest_main");
push(@cmd,"--gtest_list_tests")           if  $OPT{list};
push(@cmd,"--gtest_filter=$OPT{filter}")  if  $OPT{filter};
push(@cmd,"--gtest_repeat=$OPT{repeat}")  if  $OPT{repeat};
push(@cmd,"--gtest_shuffle")              if  $OPT{shuffle};

exec(@cmd);
