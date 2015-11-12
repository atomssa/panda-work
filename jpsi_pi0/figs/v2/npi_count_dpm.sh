
echo "{" > count.12.C
echo "TH1F* h=new TH1F(\"h\",\"h\",11,0.5,11.5);" >> count.12.C
for i in {3..10}; do echo -n "h->Fill("$i", " >> count.12.C; c=$(grep ^#$i count.12.txt | awk 'BEGIN{s=0}{s+=$5}END{print s}'; ); echo $c");" >> count.12.C; done
echo "h->Draw();" >> count.12.C
echo "}" >> count.12.C
root count.12.C

echo "{" > count.5.C
echo "TH1F* h=new TH1F(\"h\",\"h\",11,0.5,11.5);" >> count.5.C
for i in {3..10}; do echo -n "h->Fill("$i", " >> count.5.C; c=$(grep ^#$i count.5.txt | awk 'BEGIN{s=0}{s+=$5}END{print s}'; ); echo $c");" >> count.5.C; done
echo "h->Draw();" >> count.5.C
echo "}" >> count.5.C
root count.5.C

echo "{" > count.8.C
echo "TH1F* h=new TH1F(\"h\",\"h\",11,0.5,11.5);" >> count.8.C
for i in {3..10}; do echo -n "h->Fill("$i", " >> count.8.C; c=$(grep ^#$i count.8.txt | awk 'BEGIN{s=0}{s+=$5}END{print s}'; ); echo $c");" >> count.8.C; done
echo "h->Draw();" >> count.8.C
echo "}" >> count.8.C
root count.8.C
