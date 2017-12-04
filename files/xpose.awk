{
  for (i=1;i<=NF;i++) a[NR,i]=$i
  if(big <= NF) big=NF
}
END {
   for(i=1;i<=big;i++)
   {
     for(j=1;j<=NR;j++) printf OFS a[j,i]; printf "\n"
   }
}
