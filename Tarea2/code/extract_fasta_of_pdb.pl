use List::MoreUtils qw(uniq);

my (@residues,$PDB_file,$name);

if(!$ARGV[0] or !$ARGV[1]){ die "# usage: $0 <PDB file>\n"; }
else{ $PDBfile = $ARGV[0]; $name = $ARGV[1]; }

my %triple_to_one = (
    'CYS' , 'C', 'ASP' , 'D', 'SER' , 'S', 'GLN' , 'Q', 'LYS' , 'K',
    'ILE' , 'I', 'PRO' , 'P', 'THR' , 'T', 'PHE' , 'F', 'ASN' , 'N', 
    'GLY' , 'G', 'HIS' , 'H', 'LEU' , 'L', 'ARG' , 'R', 'TRP' , 'W', 
    'ALA' , 'A', 'VAL' , 'V', 'GLU' , 'E', 'TYR' , 'Y', 'MET' , 'M'
);

# Leer Residuos de aminoácidos.
open(PDB, $PDBfile) || die "# $0 : No puedo leer $PDBfile\n";
while(<PDB>)
{
	last if(/^ENDMDL/); # para estructuras NMR como 1lfu, TER es otra opcion
	next if not(/^ATOM\s+\d+\s+\w+\s+(\w{3})\s+\w+\s+(\d+)/); # Extracion de Residuo y su numero de sequencia.
	push(@residues, "$1\t$2");
}
close(PDB);

@residues = map {$triple_to_one{substr($_, 0 ,3)}} uniq @residues; # Convert each residue to single letter code.
print ">" . $name . "\n" .join("", @residues) . "\n"; # Print as fasta format.

