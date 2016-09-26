awk 'BEGIN{Ma=0; Mb=0; Sa=0; Sb=0; Covab=0;}
                        {a[NR]=$1;  b[NR]=$2; c[NR]=($1-$2)*($1-$2); label[NR]=$3;
                         Ma=Ma+a[NR];  Mb=Mb+b[NR];}
END{ 
eg=0; for(i=1;i<=NR-1;i++){if(a[i]!=a[i+1]){eg++}}; 
if(eg!=0){
if(NR>=2){for(i=1;i<=NR;i++){Sa=Sa+(a[i]-Ma/NR)*(a[i]-Ma/NR)/NR;
                        Sb=Sb+(b[i]-Mb/NR)*(b[i]-Mb/NR)/NR;
                        e=e+c[i]/NR;
                        Covab=Covab+a[i]*b[i]-Ma/NR*Mb/NR;}
for(i=1;i<=NR;i++){     Za[i]=(a[i]-Ma/NR)/sqrt(Sa);
                        Zb[i]=(b[i]-Mb/NR)/sqrt(Sb);}
for(i=1;i<=NR;i++){     C=C+Za[i]*Zb[i]/NR; }
                        B=(Covab/NR)/(Sa);
                        A=Mb/NR-B*Ma/NR;
error=sqrt(e);
print "#", "'$1'",   C, error, A, B, NR;
{ for(i=1;i<=NR;i++){ print a[i], (b[i]-A)/B,label[i] }};}}
if(eg==0){print "#", "'$1'",  0, "N", "N", "N", NR;}
}' $1
