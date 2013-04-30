
BEGIN { ex=90; }
{t1=$4-ex; t2=$5+ex-1; if(t1<0) {t1 = 1; }  printf "%s\t.\t.\t%d\t%d\t%s\t.\t.\t%s\n",  $1,t1,t2,$6,$9 }
