#WBL 2 August 2022 relable GP random_tree() in bash expr format

(FNR==1){
  print "#text_expr.awk $Revision: 1.5 $", FILENAME;
}
(index($0,"Revision")){print "#"$0}
($1=="seed"){srand($2);print "#"$1,$2}
#(index($0,"%")==1||$1=="evalstacksize:"||$1=="seed"){print "#"$0;}
#($1=="Created"){printf("%s %s %s %s ",$5,$6,$9,$10);}
#($1=="ran" || $2=="took"){print "#"$0}

($1=="avx" && (s=index($0,"="))){
  printf("\nexpr ");
  eval(2);
  printf "\n";
  next;
}
#(0 && index($0,"MUL")){
#  print $0;
#  eval(0);
#  printf "\n";
#}

#we start with just binary function MUL, but replace it with 14
function rand_function(text,  n,t) {
  ;
  n = split("+ - * / % | & < <= = == != >= >",t);
  return t[1+int(n * rand())];
}
function Print(text) {
  ;
  if(index(text,"MUL")) printf("\"%s\" ",rand_function())
  else                  printf("%d ",int(1+32768*rand()));
}
function eval(I,  I_,nargs,i,simple){
  I++;
  I_ = I;
  nargs = index($I,"MUL")? 2 : 0;
  if(nargs==2) {
    simple = index($(I+1),"MUL") == 0;
    if(!simple) printf("\"(\" ");
    I = eval(I);
    if(!simple) printf("\")\" ");
    Print($I_);
    simple = index($(I+1),"MUL") == 0;
    if(!simple) printf("\"(\" ");
    I = eval(I);
    if(!simple) printf("\")\" ");
  } else Print($I);
  return I;
}
