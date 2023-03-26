OUTDIR=analysis

## CIRCOS
# extract records and protein sequences from plasmids' GenBank files
python scripts/gbkNPextractor.py -i plasmids -o $OUTDIR/plasmids.np -p -r -c -q record.id,protein_id,product -u '_' -x '|'
# combine proteins from plasmids into a single file
cat $OUTDIR/plasmids.np/prot/* > $OUTDIR/proteins.fasta
# combine nucleotide sequences of plasmids into a single file
cat $OUTDIR/plasmids.np/records/* > $OUTDIR/records.fasta
# generate Circos plots based on both nucleoite sequences
python scripts/circos_wrapper.py --nucl $OUTDIR/records.fasta --name  ANT_H3 --output $OUTDIR/circos --order plasmids.nucl.labels.tsv --add_nucl_label plasmids.nucl.labels.tsv --i_label_size 10 --image_size 1000

## BLASTP
# run blastp on proteins from plasmids
PROTBOUT=$OUTDIR/ANT_H3.blastp.e1e-10_p30_q90.tsv
blastp -query $OUTDIR/proteins.fasta -subject $OUTDIR/proteins.fasta -evalue 1e-10 -qcov_hsp_perc 90 -out $PROTBOUT -outfmt "6 std qcovhsp"

## NETWORK
# make a directory for network files
mkdir $OUTDIR/network
# create edges and nodes files from blastp output
python scripts/bout2ENdict.py --infile $PROTBOUT --outfile $OUTDIR/network/net_noG --preset qcovhsp -fhec -id --addpq --att 'Node,Acc,Prot_ID,Product,Coords' --grouping 0 -p 30
# match plasmids' descriptions with nodes
python scripts/column_matcher.py -r $OUTDIR/network/net_noG.nodes -s plasmids.desc.tsv -o $OUTDIR/network/net_noG.nodes.upt -c 2 -C 2 --header_prefix '#'
# create graphml network file from edges and nodes files
python scripts/anNe.py --nodes $OUTDIR/network/net_noG.nodes.upt --edges $OUTDIR/network/net_noG.edges --outfile $OUTDIR/network/net_noG