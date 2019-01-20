use warnings;
use strict;
use Getopt::Long;
use Cwd;


sub parse_metadata {
    my $opts = shift;
    my %run = ();
    open my $fh, "<", $opts->{"input"} || die $!;
    while(<$fh>){
        chomp;
        next if $_=~m/^\s*$/;
        if($opts->{"paired"}){
            my ($sample, $read1, $read2) = split /\t/;
            $run{$sample} = [$read1, $read2];
        }else{
            my ($sample, $read) = split /\t/;
            $run{$sample} = [$read]; 
        }
    }
    close $fh;
    return \%run;
}

sub parse_config {
    
    my $opts = shift;
    my %conf = ();
    
    open my $fh, "<", $opts->{"config"} || die $!;
    while(<$fh>){
        chomp;
        next if $_ =~ m/^\s*$/;
        my ($key, $val) = split /\t/;
        $conf{$key} = $val;
    }
    close $fh;
    return \%conf;
}

sub show_help {
    print "\n  --input <input metadata>\n",
    "  --output <output directory>\n",
    "  --thread <number of threads>\n",
    "  --log <keep log>\n",
    "  --verbose <verbose>\n",
    "  --config <configuration file>\n",
    "  --paired <paired end reads>\n",
    "  --help <show help>\n";
    return 0;
}

sub check_input {
    my $opts = shift;
    if($opts->{"help"}){
        show_help;
        return 0;
    }
    if(!defined $opts->{"input"} || !defined $opts->{"output"} || !defined $opts->{"config"}){
        die "Input/Output/config options not valid!\n";
    }
    unless(-r $opts->{"input"} && -r $opts->{"config"}){
        die "Input/config file not readable!\n";
    }
    return 0;
}


sub get_input {

    my %opts = ("input" => undef,
                "output" => undef,
                "config" => undef,
                "thread" => 1,
                "paired" => 0,
                "log" => 0,
                "verbose" => 0,
                "help" => 0);

    GetOptions("verbose" => \$opts{"verbose"},
                "log" => \$opts{"log"},
                "thread=i" => \$opts{"thread"}, 
                "input=s" => \$opts{"input"}, 
                "output=s" => \$opts{"output"},
                "config=s" => \$opts{"config"},
                "paired" => \$opts{"paired"},
                "help" => \$opts{"help"});
    check_input(\%opts);
    return \%opts;
}





## Main ##
my ($opts, $conf, $run);
eval {$opts = get_input();};
if($@){
    warn $@;
    exit 1;
}
if($opts->{"help"} == 0){
    print "*"x20, "\n";
    map{printf "* %-10s:%s\n", $_, $opts->{$_}}(keys %{$opts});
    $conf = parse_config($opts);
    map{printf "* %-10s:%s\n", $_, $conf->{$_}}(keys %{$conf});
    $run = parse_metadata($opts);
    map{print "* ", $_,"\t",join(",", @{$run->{$_}}), "\n"}(keys %{$run});
    print "*"x20, "\n";
}else{
    exit 1;
}


print "\nChecking Clark Files..\n";
if(!-e $conf->{"CLARK"} || !-e $conf->{"DB"}){
    print "Clark executable or database not found\n";
    print "Prepare CLARK Environment ? (y/N): ";
    my $req = <>;
    if($req=~m/y/i){
        if(system("wget -O CLARK.tar.gz http://clark.cs.ucr.edu/Download/CLARKV1.2.4.tar.gz")){
            print "Failed to download clark source\n";
            exit 0;
        }else{
            if(system("tar -xvzf CLARK.tar.gz")){
                print "Failed to open clark archive\n";
                exit 0;
            }else{
                my @clark_dir_name = grep{-d $_ && $_=~m/clarks.+/i} glob("*");
                my $clark_dir_name = $clark_dir_name[0];
                print "Enter directory path to download clark db: ";
                chomp(my $clark_dir = <>);
                printf("--> bash %s/set_targets.sh %s bacteria viruses\n", $clark_dir_name, $clark_dir);
                if(system(sprintf("bash %s/set_targets.sh %s bacteria viruses", $clark_dir_name, $clark_dir))){
                    print "Failed to initialize clark\n";
                }else{
                    $conf->{"CLARK"} = $clark_dir_name;
                    $conf->{"DB"} = $clark_dir;
                }
            }
        }
    }else{
        exit 1;
    }
}

print "Preparing Files for Clark..\n";
if($opts->{"paired"}){
    open OUTF, ">ClarkF.txt" || die $!; 
    open OUTR, ">ClarkR.txt" || die $!;
    open RES, ">ClarkRes.txt" || die $!;

    if(chmod("0644", "ClarkF.txt", "ClarkR.txt", "ClarkRes.txt") != 3){
        print "Failed to set file permissions\n";
    }
    print "Writing reads ClarkF.txt,ClarkR.txt,ClarkRes.txt\n";
    foreach my $sample (keys %{$run}){
        print OUTF $run->{$sample}->[0], "\n";
        print OUTR $run->{$sample}->[1], "\n";
        print RES $opts->{"output"}.'/'.$sample, "\n";
    }
    close OUTF;
    close OUTR;
    close RES;

    my $cur_dir = cwd();
    chdir $conf->{"CLARK"};
    print "Running Clark..\n";
    printf("--> %s -m 0 -n %d --light -P %s %s -R %s\n", $conf->{"CLARK"}.'/classify_metagenome.sh', $opts->{"thread"}, $cur_dir."/ClarkF.txt", $cur_dir."/ClarkR.txt", $cur_dir."/ClarkRes.txt");
    system(sprintf("bash %s -m 0 -n %d --light -P %s %s -R %s", $conf->{"CLARK"}.'/classify_metagenome.sh', $opts->{"thread"}, $cur_dir."/ClarkF.txt", $cur_dir."/ClarkR.txt", $cur_dir."/ClarkRes.txt"));
}else{
    open OUTF, ">ClarkF.txt" || die $!; 
    open RES, ">ClarkRes.txt" || die $!;
    if(chmod("0644", "ClarkF.txt", "ClarkRes.txt") != 2){
        print "Failed to set file permissions\n";
    }
    print "Writing reads ClarkF.txt,ClarkRes.txt\n";
    foreach my $sample (keys %{$run}){
        print OUTF $run->{$sample}->[0], "\n";
        print RES $opts->{"output"}.'/'.$sample, "\n";
    }
    close OUTF;
    close RES;    

    my $cur_dir = cwd();
    chdir $conf->{"CLARK"};
    print "Running Clark..\n";
    printf("--> %s -m 0 -n %d --light -O %s -R %s\n", $conf->{"CLARK"}.'/classify_metagenome.sh', $opts->{"thread"}, $cur_dir."/ClarkF.txt", $cur_dir."/ClarkRes.txt");
    system(sprintf("bash %s -m 0 -n %d --light -O %s -R %s", $conf->{"CLARK"}.'/classify_metagenome.sh', $opts->{"thread"}, $cur_dir."/ClarkF.txt", $cur_dir."/ClarkRes.txt"));
}

## Estimate Abundance ##
foreach my $sample (keys %{$run}){
    print "Estimating Abundance for sample $sample\n";
    printf("--> %s --highconfidence -D %s -F %s > %s\n", $conf->{"CLARK"}.'/estimate_abundance.sh', $conf->{"DB"}, $opts->{"output"}.'/'.$sample.'.csv', $opts->{"output"}.'/Abundance_'.$sample.'.csv');
    system(sprintf("%s --highconfidence -D %s -F %s > %s\n", $conf->{"CLARK"}.'/estimate_abundance.sh', $conf->{"DB"}, $opts->{"output"}.'/'.$sample.'.csv', $opts->{"output"}.'/Abundance_'.$sample.'.csv'))
}





exit 0;