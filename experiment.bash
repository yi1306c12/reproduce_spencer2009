STR=$1
DIRNAME='./single/'
LOGNAME=$DIRNAME$STR'_log.txt'

GDL_STARTUP=setup.txt gdl -quiet -e "ctx, '$STR', "$STR > $LOGNAME