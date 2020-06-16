csi_image_v6   << EOF
$1.cmplx
w
dummy_mrs_mask.real
0.8
t
$1b
$2.peak
1
21
0 5
n
n
n
n
s
$1b_cor.cmplx
0 0 2 4
EOF

cp $1.ddf $1b_cor.ddf
cp $1.ddf $1b_cor_back.ddf
cp $1.ddf $1b_cor_sum.ddf
mkdir $1b_cor_met
mv $1b*0*.real $1b_cor_met/.
mv $1b*0*.idf $1b_cor_met/.
mv $1b*tfr.real $1b_cor_met/.
mv $1b*tfr.idf $1b_cor_met/.
mv $1b*tph.real $1b_cor_met/.
mv $1b*tph.idf $1b_cor_met/.
mv $1b*tnsd.real $1b_cor_met/.
mv $1b*tnsd.idf $1b_cor_met/.
mv $1b*tpar.real $1b_cor_met/.
mv $1b*tpar.idf $1b_cor_met/.
mv $1b*trefh.idf $1b_cor_met/.
mv $1b*trefh.real $1b_cor_met/.
mv $1b*wn*.idf $1b_cor_met/.
mv $1b*wn*.real $1b_cor_met/.
mv $1b*.tab $1b_cor_met/.

