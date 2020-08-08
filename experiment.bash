VARNAME=$1
RATE=$2
STR=$VARNAME$RATE
DIRNAME='./single/'
LOGNAME=$DIRNAME'log_'$STR'.txt'

EXPRESSION='ctx, "'$STR'", '$VARNAME'='$RATE
echo $EXPRESSION > $LOGNAME
echo $PATH
GDL_STARTUP=setup.txt gdl -e $EXPRESSION > $LOGNAME
#GDL_STARTUP=setup.txt gdl -e 'ctx, "''", con_t='$RATE > ./single/con_t$RATE.log