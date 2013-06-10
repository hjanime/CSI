
BEGIN { ex=90; }
{t1=$4-ex; t2=$5+ex-1; if(t1<0) {t1 = 1; }  printf "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n",  $1,$2,$3,t1,t2,$6,$7,$8,$9 }
