for j in `cat $1  | awk '{print $1}' | uniq`; 
do 
#awk '($1=="'$j'"){print "'$j'", $2}' $1 | sort -n -k 2 | tail -1 | sed 's/-/ -/g' | awk '{print $1+(-$1-$2)/2,$3}'
awk '($1=="'$j'"){print "'$j'", $2}' $1 | sort -n -k 2 | tail -1 | sed 's/-/ -/g' | awk '{print $1,$3}'

#; 
done  | awk '{s=s+$2; a[NR]=$1; b[NR]=$2} END{for(i=1;i<=NR;i++){f=f+(b[i]-s/NR)^2} for(i=1;i<=NR;i++){print a[i], (b[i]-s/NR)/sqrt(f/NR)}}' | awk '{a[NR]=$1; z[NR]=$2}  END{for(i=1;i<=NR;i++){
if(i==1){print a[i],  z[i]/3}
if(i==2){print a[i], (z[i]+z[i-1]+z[i+1])/5}
if(i==3){print a[i], (z[i]+z[i-1]+z[i+1]+z[i+2]+z[i-2])/7}
if(i==NR)  {print a[i], z[i]/3}
if(i==NR-1){print a[i], (z[i]+z[i-1]+z[i+1])/5}
if(i==NR-2){print a[i], (z[i]+z[i-1]+z[i+1]+z[i+2]+z[i-2])/7}

if((i>3)&&(i<NR-2)){print (a[i-3]+a[i-2]+a[i-1]+a[i]+a[i+1]+a[i+2]+a[i+3])/7, (z[i-3]+z[i-2]+z[i-1]+z[i]+z[i+1]+z[i+2]+z[i+3])/7,0} }}' | awk '{s=s+$2; a[NR]=$1; b[NR]=$2} END{for(i=1;i<=NR;i++){f=f+(b[i]-s/NR)^2} for(i=1;i<=NR;i++){print a[i], (b[i]-s/NR)/sqrt(f/NR)}}'
# | awk '($2>0)' 

