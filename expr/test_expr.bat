#WBL 2 August 2022 test coreutils 9.1 expr

#seed     $1

set seed=171534
if( $1 != "") set seed=$1

rm -rf /tmp/test_expr
mkdir -p /tmp/test_expr
if($status) exit $status;

echo $0 '$Revision: 1.11 $'  $seed `pwd` `date` $HOST

set ok=0
set len=1
while ( $len < 400 )
  set i=1
  while ( $i <= 60 )
    set ll=(`printf "%03d" $len`)
    set ii=(`printf "%02d" $i`)
    ./bench_rand         pMaxExpr:$len pPopSeed:$seed | gawk -f test_expr.awk > ! /tmp/test_expr/test{$ll}_{$ii}.bat
    setenv save $status;
    if($save) then
      echo "./bench_rand pMaxExpr:$len pPopSeed:$seed | gawk -f test_expr.awk, status $save";
      exit $save;
    endif
    chmod +x /tmp/test_expr/test{$ll}_{$ii}.bat
    /tmp/test_expr/test{$ll}_{$ii}.bat
    if($status) then
      rm -f /tmp/test_expr/test{$ll}_{$ii}.bat
    else
      set ok=(`expr $ok + 1`)
    endif
    set i=(`expr $i + 1`)
    set seed=(`expr $seed + 99`)
  end
  echo "$0 done $len, $ok so far"
  set len=(`expr $len + 2`)
end

echo "$0 done $ok $len $seed" `date`

