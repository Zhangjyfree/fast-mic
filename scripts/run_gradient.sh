#!/usr/bin/env bash
# 跑 fast-mic 穿过 8 级培养基梯度（akk_minimal → western_diet_mucin），
# 收集每级的 pairwise 互作分布，供下游 analyze_gradient.py 汇总。
#
# Run fast-mic across the 8-level medium gradient (akk_minimal →
# western_diet_mucin) and collect pairwise interaction-type distributions
# per level for downstream aggregation.
#
# 用法 / Usage:
#   bash scripts/run_gradient.sh \
#        --group1 test/akk/akk_gapseq_xml \
#        --group2 test/UHGG_v2/final_bacteria_gapseq_xml \
#        --threads 12 --full-tsv
#
#   bash scripts/run_gradient.sh \
#        --group1 test/lac/lac_gapseq_xml \
#        --group2 test/UHGG_v2/final_bacteria_gapseq_xml \
#        --threads 12 --full-tsv \
#        --out results/lac_vs_uhgg          # ← 显式指定输出目录
#
# 选项 / Options:
#   --models DIR        单组 all-vs-all
#   --group1 DIR        跨组第一组
#   --group2 DIR        跨组第二组
#   --threads N         线程数（0 = 全核）
#   --full-tsv          额外输出 *.full.tsv（逐代谢物交叉营养）
#   --out DIR           输出目录（默认从 group 名自动推导，见下）
#
# 输出目录 / Output:
#   若不指定 --out，默认自动推导，避免不同分析互相覆盖：
#     group1=akk, group2=UHGG_v2  →  gradient_result_akk_vs_UHGG_v2/
#     group1=lac, group2=UHGG_v2  →  gradient_result_lac_vs_UHGG_v2/
#     --models akk                →  gradient_result_akk/
#   目录内容:
#     L0_base.tsv ... L9_mos.tsv
#     L*.full.tsv   (若启用 --full-tsv)
#     gradient_run.log

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")"/.. && pwd)"
BIN="$SCRIPT_DIR/target/release/fast-mic"
MEDIA_DIR="$SCRIPT_DIR/media"
OUT_DIR=""                       # 空 = 自动推导（见下方 derive 逻辑）

# ── 命令行选项 / CLI options ─────────────────────────────────────────────────
MODELS=""
GROUP1=""
GROUP2=""
THREADS=0
WANT_FULL_TSV=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --models)      MODELS="$2"; shift 2 ;;
        --group1)      GROUP1="$2"; shift 2 ;;
        --group2)      GROUP2="$2"; shift 2 ;;
        --threads)     THREADS="$2"; shift 2 ;;
        --out)         OUT_DIR="$2"; shift 2 ;;
        --full-tsv)    WANT_FULL_TSV=1; shift ;;
        -h|--help)
            sed -n '2,40p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ── 从模型目录路径推导一个简短标签 / derive a short label from a model dir ──────
# 如果叶子目录名很泛（含 gapseq / xml / bacteria / final），用其父目录名。
#   test/akk/akk_gapseq_xml            -> akk
#   test/UHGG_v2/final_bacteria_..._xml -> UHGG_v2
#   test/lac/lac_gapseq_xml            -> lac
group_label() {
    local d leaf parent
    d="${1%/}"                         # 去掉结尾斜杠
    leaf="$(basename "$d")"
    parent="$(basename "$(dirname "$d")")"
    if [[ "$leaf" =~ (gapseq|xml|bacteria|final|sbml) ]]; then
        echo "$parent"
    else
        echo "$leaf"
    fi
}

# ── 自动推导 OUT_DIR（仅当用户未用 --out 显式指定时）/ auto-derive OUT_DIR ──────
if [[ -z "$OUT_DIR" ]]; then
    if [[ -n "$MODELS" ]]; then
        OUT_DIR="$SCRIPT_DIR/gradient_result_$(group_label "$MODELS")"
    elif [[ -n "$GROUP1" && -n "$GROUP2" ]]; then
        OUT_DIR="$SCRIPT_DIR/gradient_result_$(group_label "$GROUP1")_vs_$(group_label "$GROUP2")"
    else
        OUT_DIR="$SCRIPT_DIR/gradient_result"   # fallback（校验阶段会报错）
    fi
fi

# 校验 / Validate
if [[ -z "$MODELS" && ( -z "$GROUP1" || -z "$GROUP2" ) ]]; then
    echo "Error: provide either --models DIR or both --group1 DIR --group2 DIR" >&2
    exit 1
fi
if [[ ! -x "$BIN" ]]; then
    echo "Error: fast-mic binary not found at $BIN — run 'cargo build --release' first." >&2
    exit 1
fi
if ! ls "$MEDIA_DIR"/gradient_L*_gapseq.csv >/dev/null 2>&1; then
    echo "Error: gradient media not found in $MEDIA_DIR. Run scripts/generate_medium_gradient.py first." >&2
    exit 1
fi

mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/gradient_run.log"
: > "$LOG"
echo "===========================================================" | tee -a "$LOG"
echo " Gradient run started $(date -u +%FT%TZ)" | tee -a "$LOG"
echo " OUTPUT DIR ▶  $OUT_DIR" | tee -a "$LOG"
echo "===========================================================" | tee -a "$LOG"

# ── 物种集摘要 / Species-set summary ─────────────────────────────────────────
if [[ -n "$MODELS" ]]; then
    N=$(find "$MODELS" -maxdepth 1 \( -name '*.xml' -o -name '*.sbml' \) | wc -l | tr -d ' ')
    echo "Mode: all-vs-all on $MODELS  ($N models, $((N*(N-1)/2)) pairs per medium)" | tee -a "$LOG"
else
    N1=$(find "$GROUP1" -maxdepth 1 \( -name '*.xml' -o -name '*.sbml' \) | wc -l | tr -d ' ')
    N2=$(find "$GROUP2" -maxdepth 1 \( -name '*.xml' -o -name '*.sbml' \) | wc -l | tr -d ' ')
    echo "Mode: cross-group $GROUP1 x $GROUP2  ($N1 x $N2 = $((N1*N2)) pairs per medium)" | tee -a "$LOG"
fi
echo "Threads: $THREADS" | tee -a "$LOG"
echo "Output:  $OUT_DIR" | tee -a "$LOG"
echo "" | tee -a "$LOG"

# ── 主循环 / Main loop over the 8 gradient media ────────────────────────────
for med in "$MEDIA_DIR"/gradient_L*_gapseq.csv; do
    # 提取标签 / extract short label
    label=$(basename "$med" .csv | sed 's/^gradient_//;s/_gapseq$//')
    out="$OUT_DIR/${label}.tsv"
    full="$OUT_DIR/${label}.full.tsv"

    n_cpd=$(($(wc -l < "$med") - 1))
    echo "[$(date +%H:%M:%S)] Running $label ($n_cpd compounds) → $out" | tee -a "$LOG"

    args=( --medium-file "$med" --threads "$THREADS" -o "$out" )
    if [[ $WANT_FULL_TSV -eq 1 ]]; then
        args+=( --full-tsv "$full" )
    fi
    if [[ -n "$MODELS" ]]; then
        # 让 fast-mic 自动遍历目录里的所有 SBML / let fast-mic glob the directory
        # shellcheck disable=SC2206
        models_arr=( "$MODELS"/*.xml "$MODELS"/*.sbml )
        # 过滤掉 glob 不到的 fallback / drop unmatched globs
        models_clean=()
        for m in "${models_arr[@]}"; do [[ -f "$m" ]] && models_clean+=( "$m" ); done
        args+=( --summary "${models_clean[@]}" )
    else
        args+=( --group1 "$GROUP1" --group2 "$GROUP2" --summary )
    fi

    # 记录到日志的同时计时 / time + log
    /usr/bin/time -p "$BIN" "${args[@]}" >>"$LOG" 2>&1
    echo "    done." | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "=== Gradient run finished $(date -u +%FT%TZ) ===" | tee -a "$LOG"
echo ""
echo "Next step:"
echo "  python3 scripts/analyze_gradient.py --results $OUT_DIR"
