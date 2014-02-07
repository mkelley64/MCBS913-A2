use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'FastaORFUtils',
    license             => 'perl',
    dist_abstract       => 'Bioinformatics lab utilities',
    dist_author         => 'Mark Kelley <meb223@unh.edu>',
    build_requires => {
        'Test::More' => '0.10',
    },
);

$builder->create_build_script();
