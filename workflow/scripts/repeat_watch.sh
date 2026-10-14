#!/bin/bash
# repeat_watch.sh -- draws the repeat figures (repeat_context.py) for each sample as soon as reproducible_phd
# publishes that species' CARP annotation (refs/{species}_repeats/MANIFEST.txt), and again if it is re-published.
# Start it in tmux on the login node:
#   tmux new-session -d -s repeat_watch 'bash workflow/scripts/repeat_watch.sh; exec bash'
# Looks every 10 min; stops when every sample has figures newer than its annotation (or a failure newer than
# both the annotation and the script), or after 3 days. Its log: logs/repeat_watch.log
CO=/netscratch/dep_mercier/grp_marques/Aaryan/CO_smk
REFS=/netscratch/dep_mercier/grp_marques/Aaryan/refs
PY=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/python
SAMPLES="Dparadoxa_split:Drosera_paradoxa Dbinata_hap1:Drosera_binata"
EVERY=${REPEAT_WATCH_EVERY:-600}
SCRIPT=workflow/scripts/repeat_context.py
cd "$CO" || exit 1
mkdir -p logs reports/repeat_context
end=$(( $(date +%s) + 3 * 86400 ))
echo "$(date '+%F %T') watching for: $SAMPLES (every $EVERY s)" | tee -a logs/repeat_watch.log
while [ "$(date +%s)" -lt "$end" ]; do
    left=0
    for x in $SAMPLES; do
        s=${x%%:*}
        sp=${x##*:}
        man=$REFS/${sp}_repeats/MANIFEST.txt
        png=reports/repeat_context/repeat_landscape_$s.png
        bad=logs/repeat_context_$s.failed
        if [ ! -s "$man" ]; then
            left=$((left + 1))
            continue
        fi
        [ -s "$png" ] && [ "$png" -nt "$man" ] && continue
        [ -e "$bad" ] && [ "$bad" -nt "$man" ] && [ "$bad" -nt "$SCRIPT" ] && continue
        echo "$(date '+%F %T') $sp annotation published: repeat_context.py $s" | tee -a logs/repeat_watch.log
        if "$PY" "$SCRIPT" "$s" > "logs/repeat_context_$s.log" 2>&1; then
            rm -f "$bad"
            grep '^wrote' "logs/repeat_context_$s.log" | sed 's/^/   /' | tee -a logs/repeat_watch.log
        else
            touch "$bad"
            echo "   FAILED, last lines of logs/repeat_context_$s.log:" | tee -a logs/repeat_watch.log
            tail -n 6 "logs/repeat_context_$s.log" | sed 's/^/   /' | tee -a logs/repeat_watch.log
        fi
    done
    if [ $left -eq 0 ]; then
        echo "$(date '+%F %T') every sample is drawn; the watcher stops" | tee -a logs/repeat_watch.log
        exit 0
    fi
    sleep "$EVERY"
done
echo "$(date '+%F %T') gave up after 3 days" | tee -a logs/repeat_watch.log
