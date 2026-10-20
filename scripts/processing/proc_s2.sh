#!/usr/bin/env bash
# S2 robustness: full download, random 4M (same depth as first-4M), into quant_s2/
# usage: proc_fast.sh RUN READCOUNT "url1;url2" "md5_1;md5_2" INDEX BYTES
# Small files (<1.2 GB total): full MD5-checked download + random subsample of 8M reads (fixed seed).
# Large files: first 4M reads / read pairs streamed (both mates read in lock-step; no partial restarts).
export PATH=${SALMON_BIN:-/opt/salmon-1.10.0/bin}:$PATH
run=$1; rc=$2; urls=$3; md5s=$4; index=$5; bytes=$6
out=${PROJECT_DIR:-.}/quant_s2/$run; tmp=${PROJECT_DIR:-.}/fq/s2_$run; L=${PROJECT_DIR:-.}/logs
[[ -f $out/quant.sf ]] && exit 0
mkdir -p "$tmp"
IFS=';' read -ra U <<< "$urls"; IFS=';' read -ra H <<< "$md5s"
total=$(python3 -c "print(sum(int(x) for x in '$bytes'.split(';') if x))")
files=()
if true; then
  mode=random8M
  for i in "${!U[@]}"; do
    f=$tmp/$(basename "${U[$i]}"); files+=("$f")
    for a in 1 2 3 4; do
      [[ -f $f.ok ]] && break
      curl -sS -C - --retry 5 --retry-delay 15 -m 7200 -o "$f" "https://${U[$i]}" 2>>$L/$run.dl
      [[ "$(md5sum "$f" | cut -d' ' -f1)" == "${H[$i]}" ]] && { touch "$f.ok"; break; }
      (( a >= 2 )) && rm -f "$f"
      sleep 15
    done
    [[ -f $f.ok ]] || { echo "$run FAILED download"; exit 1; }
  done
  p=$(python3 -c "print(min(1.0, 4000000/max(1,$rc)))")
  need=$(python3 -c "print(int(min(4000000,$rc)*0.95))")
else
  mode=first4M
  need=$(python3 -c "print(int(min(4000000,$rc)*0.99))")
  IFS=';' read -ra B <<< "$bytes"
  nl=$(python3 -c "print(4*min(4000000,$rc))")
  for i in "${!U[@]}"; do
    f=$tmp/part_$i.fastq.gz; files+=("$f")
    for factor in 1.12 1.35 1.7; do
      want=$(python3 -c "print(min(int(${B[$i]})-1, int(int(${B[$i]})*4000000*$factor/max(1,$rc))))")
      for a in $(seq 1 15); do
        have=$(stat -c %s "$f" 2>/dev/null || echo 0)
        (( have > want )) && break
        curl -fsS -m 1800 -r ${have}-${want} "https://${U[$i]}" >> "$f" 2>>$L/$run.dl || sleep 10
      done
      got=$(gzip -dc "$f" 2>/dev/null | head -n $nl | wc -l)
      (( got >= nl )) && { touch "$f.ok"; break; }
      echo "$run part $i short ($got lines) at factor $factor"
      mv -f "$f" "$f.bad"   # corrupted or too short: start this part fresh with a larger margin
    done
    [[ -f $f.ok ]] || { echo "$run FAILED partial download"; exit 1; }
  done
fi
# deterministic per-record sampling (identical selection in both mate files; mawk srand() is not reproducible across processes)
src_rand() { gzip -dc "$1" | mawk -v p="$p" 'NR%4==1{i=(NR-1)/4; x=i*0.6180339887498949+0.5; k=((x-int(x))<p)} k'; }
src_head() { gzip -dc "$1" 2>/dev/null | head -n 16000000; }
for attempt in 1 2 3; do
  rm -rf "$out.tmp"
  if [[ $mode == random8M ]]; then
    if [[ ${#files[@]} -ge 2 ]]; then
      salmon quant -i "$index" -l A -1 <(src_rand "${files[0]}") -2 <(src_rand "${files[1]}") -p 2 \
        --validateMappings --seqBias --gcBias -o "$out.tmp" -q 2> $L/$run.err
    else
      salmon quant -i "$index" -l A -r <(src_rand "${files[0]}") -p 2 --validateMappings --seqBias \
        -o "$out.tmp" -q 2> $L/$run.err
    fi
  else
    if [[ ${#files[@]} -ge 2 ]]; then
      salmon quant -i "$index" -l A -1 <(src_head "${files[0]}") -2 <(src_head "${files[1]}") -p 2 \
        --validateMappings --seqBias --gcBias -o "$out.tmp" -q 2> $L/$run.err
    else
      salmon quant -i "$index" -l A -r <(src_head "${files[0]}") -p 2 --validateMappings --seqBias \
        -o "$out.tmp" -q 2> $L/$run.err
    fi
  fi
  n=$(python3 -c "import json;print(json.load(open('$out.tmp/aux_info/meta_info.json'))['num_processed'])" 2>/dev/null || echo 0)
  if (( n >= need )); then
    mv "$out.tmp" "$out"; rm -rf "$tmp"
    echo "$run ok $mode $n reads $(python3 -c "import json;print(round(json.load(open('$out/aux_info/meta_info.json'))['percent_mapped'],1))")%"
    exit 0
  fi
  echo "$run attempt $attempt incomplete ($n < $need)"
  sleep 20
done
echo "$run FAILED"
