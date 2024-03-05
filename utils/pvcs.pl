#!/usr/bin/env perl

# ICON
#
# ------------------------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#
# Luis Kornblueh, MPIM, 2012-05-31
# Luis Kornblueh, MPIM, 2013-10-16
#                 - add functions to return strings with
#                   the revision information
# Luis Kornblueh, MPIM, 2024-01-15
#                 - remove svn support
#                 - add VERSION file replacing querying
#                   git when git is not available
# Sergey Kosukhin, MPIM, 2024-03-01
#                 - replace VERSION file support with config.status information
#                 - polish messages and the generated code
#
# ------------------------------------------------------------------------------

=pod

=head1 NAME

pvcs - ICON version provenance collection tool (former Provenance collection for VCS)

=head1 DESCRIPTION

Generates the version.c source file containing information on the compiled
version of ICON (repository url, branch, tag, revision).

=head1 SYNOPSIS

pvcs [options]

=head1 OPTIONS

=over 8

=item B<--srcdir <directory>>

root source directory of ICON

=item B<--builddir <directory>>

root build directory of ICON

=item B<--help>

print the usage documentation and exit

=back

=cut

#_____________________________________________________________________________

use strict;
use warnings;

use constant { true => 1, false => 0 };

use Getopt::Long;
use Pod::Usage;

use File::Temp;
use File::Compare;
use File::Copy;

my $srcdir = '.';
my $builddir = '.';
my $help = 0;

GetOptions( 'help|?' => \$help,
            'srcdir=s' => \$srcdir,
            'builddir=s' => \$builddir,
          ) or pod2usage(-exitstatus => 2);

pod2usage(-verbose => 99, -sections => 'NAME|SYNOPSIS|OPTIONS|DESCRIPTION', -exitstatus => 1) if $help;

my $default_value = 'unknown';

my $remote_url = $default_value;
my $branch = $default_value;
my $revision = $default_value;
my $git_tag = $default_value;

my $art_branch = $default_value;
my $art_revision = $default_value;
my $art_remote_url = $default_value;

if ( -e $srcdir."/.git" ) {
    $remote_url = `git --git-dir $srcdir/.git config --get remote.origin.url`;
    chomp($remote_url);
    if (index("$remote_url", "gitlab.dkrz.de") != -1) {
        my @branches = `git --git-dir $srcdir/.git branch`;
        @branches = grep(/^\*/, @branches);
        $branch = $branches[0];
        $branch =~ s/\* *//;
        $branch =~ s/ *\n//;
        my @revisions = `git --git-dir $srcdir/.git --no-pager log --max-count=1`;
        @revisions = grep(/commit/, @revisions);
        $revision = $revisions[0];
        $revision =~ s/commit *//;
        $revision =~ s/ *\n//;
        if ( -e $srcdir."/externals/art/.git" ) {
            $art_remote_url = `git --git-dir $srcdir/externals/art/.git config --get remote.origin.url`;
            chomp($art_remote_url);
            my @art_branches = `git --git-dir $srcdir//externals/art/.git branch`;
            @art_branches = grep(/^\*/, @art_branches);
            $art_branch = $art_branches[0];
            $art_branch =~ s/\* *//;
            $art_branch =~ s/ *\n//;
            my @art_revisions = `git --git-dir $srcdir/externals/art/.git --no-pager log --max-count=1`;
            @art_revisions = grep(/commit/, @art_revisions);
            $art_revision = $art_revisions[0];
            $art_revision =~ s/commit *//;
            $art_revision =~ s/ *\n//;
        }
        $git_tag = `git describe --tags --abbrev=0 --exact-match 2>/dev/null`;
        if ($?) {
            $git_tag = $default_value;
        } else {
            $git_tag =~ s/ *\n//;
        }
    }
} elsif (-e $builddir."/config.status" ) {
    $revision = 'retrieved from the configure script';
    $git_tag = `printf '\@PACKAGE_NAME\@-\@PACKAGE_VERSION\@' | $builddir/config.status -q --file=-`;
    $git_tag =~ s/ *\n//;
} else {
    print STDERR "error: failed to collect the version information\n";
}

my $need_to_compare = false;
my $version_c;
my $fname;

if ( -e "version.c") {
    $need_to_compare = true;
    $version_c = File::Temp->new(UNLINK => false);
    $fname = $version_c->filename;
} else {
    open $version_c, ">", "version.c" or die "$0: open version.c: $!";
}

print $version_c <<"VERSION_C";
#include <string.h>

static const char remote_url[] = "$remote_url";
static const char branch[] = "$branch";
static const char revision[] = "$revision";
static const char git_tag[] = "$git_tag";
static const char art_remote_url[] = "$art_remote_url";
static const char art_branch[] = "$art_branch";
static const char art_revision[] = "$art_revision";

void repository_url(char *name, int *actual_len)
{
  if (strlen(remote_url) > *actual_len) {
    *actual_len = 0;
  } else {
    strcpy(name, remote_url);
    *actual_len = strlen(name);
  }
}

void branch_name(char *name, int *actual_len)
{
  if (strlen(branch) > *actual_len) {
    *actual_len = 0;
  } else {
    strcpy(name, branch);
    *actual_len = strlen(name);
  }
}

void revision_key(char *name, int *actual_len)
{
  if (strlen(revision) > *actual_len) {
    *actual_len = 0;
  } else {
    strcpy(name, revision);
    *actual_len = strlen(name);
  }
}

void git_tag_name(char *name, int *actual_len)
{
  if (strlen(git_tag) > *actual_len) {
    *actual_len = 0;
  } else {
    strcpy(name, git_tag);
    *actual_len = strlen(name);
  }
}

void art_repository_url(char *name, int *actual_len)
{
  if (strlen(art_remote_url) > *actual_len) {
    *actual_len = 0;
  } else {
    strcpy(name, art_remote_url);
    *actual_len = strlen(name);
  }
}

void art_branch_name(char *name, int *actual_len)
{
  if (strlen(art_branch) > *actual_len) {
    *actual_len = 0;
  } else {
    strcpy(name, art_branch);
    *actual_len = strlen(name);
  }
}

void art_revision_key(char *name, int *actual_len)
{
  if (strlen(art_revision) > *actual_len) {
    *actual_len = 0;
  } else {
    strcpy(name, art_revision);
    *actual_len = strlen(name);
  };
};
VERSION_C
close $version_c;

if ($need_to_compare) {
    if (compare($fname, "version.c") == 0) {
        unlink $fname;
    } else {
        move($fname, "version.c") or die "Copy failed: $!";
    }
}

exit 0;
