#!/bin/bash
# Nightly backup of ~/Physics/Overcharged to OneDrive (UC, remote "onedrive:") with rclone COPY (never deletes on the remote).
# Skips: deleteme* folders, __pycache__, and run directories whose solver is still running (pid file with a live process).
# Usage: backup_overcharged.sh [--dry-run]     Log: backup/logs/backup_YYYY-MM-DD.log
SRC=/home/baloewe/Physics/Overcharged; DST=onedrive:Backups/Overcharged/office; RC=/home/baloewe/bin/rclone
LOGD=$SRC/backup/logs; mkdir -p $LOGD; LOG=$LOGD/backup_$(date +%F).log
if pgrep -f "rclone copy $SRC $DST" > /dev/null; then echo "$(date) another backup is still running; skipped" >> $LOG; exit 0; fi
FILT=$(mktemp); trap 'rm -f $FILT' EXIT
{ echo "- /BACKUP_FAILED.txt"; echo "- deleteme*/**"; echo "- **/deleteme*/**"; echo "- __pycache__/**"; echo "- **/__pycache__/**"; echo "- backup/logs/**"
  for p in $(find $SRC/runs -name pid 2>/dev/null); do
      d=$(dirname $p); [ -f $d/DONE ] && continue
      pid=$(cat $p 2>/dev/null); [ -n "$pid" ] && kill -0 $pid 2>/dev/null && echo "- /${d#$SRC/}/**"
  done
  echo "+ **"; } > $FILT
echo "$(date) start ($(grep -c '^- /runs' $FILT) running run dirs skipped)" >> $LOG
$RC copy $SRC $DST --filter-from $FILT --transfers 8 --checkers 16 --onedrive-chunk-size 10M --log-file $LOG --log-level INFO --stats 10m --stats-one-line "$@"
rc=$?; echo "$(date) exit $rc" >> $LOG
# failure marker in the project root: written on any error, removed after the next clean run
MARK=$SRC/BACKUP_FAILED.txt
if [ $rc -ne 0 ]; then
    { echo "BACKUP FAILED $(date) (rclone exit $rc) — the OneDrive backup of $SRC did not complete."
      echo "Log: $LOG"
      if grep -q -i "invalid_grant\|AADSTS\|token\|unauthorized\|401" $LOG; then
          echo "Looks like an expired/revoked login. Fix: on the laptop 'ssh -N -L 53682:localhost:53682 baloewe@146.155.46.83',"
          echo "then on the office machine '~/bin/rclone config reconnect onedrive:' (browser login), then rerun $0."
      fi
      echo "Last errors:"; grep -E "ERROR|AADSTS|Failed" $LOG | tail -15; } > $MARK
elif [[ " $* " != *" --dry-run "* ]]; then
    rm -f $MARK
fi
