#!/usr/bin/awk -f

{
    if ($6 != "-") {
        $3 -=10;
        }
    else{
        $2 += 10;
        }
    printf "%s\t%d\t%d\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6;
}
