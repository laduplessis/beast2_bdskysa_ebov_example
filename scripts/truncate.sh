# Truncate tree and log files at a specified number of states
#
# Usage:
#     truncate <pattern> <nrstates>
#
# Looks for pattern*.log and pattern*.trees and truncates all files after <nrstates>
# states. Results are save in truncated/

mkdir final

# Truncate log files
for LF in `ls ${1}*.log`
do
  LINENR=`cut -f 1 $LF | grep -n ${2} | cut -f 1 -d ":"`
  head -n $LINENR $LF > truncated/${LF}
done

# Truncate tree files
for TF in `ls ${1}*.trees`
do
  LINENR=`cut -f 1 -d "=" $TF | grep -n ${2} | cut -f 1 -d ":"`
  head -n $LINENR $TF > truncated/${TF}
done
