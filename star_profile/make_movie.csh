cd movie_out
convert Power_by_Mass_*.pdf Power.jpg
mv Power-0.jpg Power-0000.jpg
mv Power-1.jpg Power-0001.jpg
mv Power-2.jpg Power-0002.jpg
mv Power-3.jpg Power-0003.jpg
mv Power-4.jpg Power-0004.jpg
mv Power-5.jpg Power-0005.jpg
mv Power-6.jpg Power-0006.jpg
mv Power-7.jpg Power-0007.jpg
mv Power-8.jpg Power-0008.jpg
mv Power-9.jpg Power-0009.jpg
mencoder mf://Power-*.jpg -mf w=1188:h=718:fps=8:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o Power.avi
rm Power-*.jpg
cd ../
