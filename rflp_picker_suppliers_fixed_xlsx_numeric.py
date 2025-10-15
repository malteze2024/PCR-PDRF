#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import sys
import re
from typing import List, Dict, Tuple, Optional
from collections import namedtuple

from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Restriction import AllEnzymes

# ====== supplier mapping (название/алиасы -> односимвольный код из Bio.Restriction) ======
_SUPPLIER_CANON = {
    "A": ["agilent"],
    "E": ["thermo fisher", "thermo scientific", "fermentas", "thermo"],
    "F": ["promega"],
    "H": ["minotech"],
    "I": ["sibenzyme", "sib enzyme"],
    "J": ["nippon gene"],
    "K": ["takara", "takara bio"],
    "M": ["roche", "roche applied science"],
    "N": ["new england biolabs", "neb"],
    "R": ["molecular biology resources", "chimerx"],
    "S": ["sigma", "sigma chemical"],
    "T": ["toyobo"],
    "V": ["vivantis"],
    "Y": ["eurx"],
    "Z": ["sinaclon", "sinaclon bioscience"],
}

def _name_to_codes(name: str) -> set[str]:
    s = (name or "").strip().lower()
    out: set[str] = set()
    if not s:
        return out
    # ровно один символ-код
    if len(s) == 1 and s.upper() in _SUPPLIER_CANON:
        return {s.upper()}
    # что-то вроде 'NEI' (без разделителей)
    if s.isalpha() and s.upper() == s and len(s) <= 8 and " " not in s:
        return {ch for ch in s if ch.upper() in _SUPPLIER_CANON}
    # алиасы/названия
    for code, aliases in _SUPPLIER_CANON.items():
        for a in aliases:
            if a in s:
                out.add(code)
                break
    return out

def parse_suppliers_arg(s: Optional[str]) -> set[str]:
    if not s:
        return set()
    tokens = [t.strip() for t in re.split(r"[,\|;]", s) if t.strip()]
    codes: set[str] = set()
    for t in tokens:
        codes |= _name_to_codes(t)
    return {c.upper() for c in codes if c.upper() in _SUPPLIER_CANON}

def enzyme_supplier_codes(e) -> set[str]:
    """Надёжно получаем коды производителей для фермента в разных версиях BioPython."""
    items = []
    # предпочтительно supplier_list()
    try:
        items = list(e.supplier_list())
    except Exception:
        pass
    # запасные варианты
    if not items:
        val = getattr(e, "suppliers", None)
        try:
            if callable(val):
                val = val()
        except Exception:
            val = None
        if isinstance(val, dict):
            items = list(val.keys()) + list(val.values())
        elif isinstance(val, (list, tuple, set)):
            items = list(val)
        elif isinstance(val, str) and val.strip():
            items = [val]

    codes: set[str] = set()
    for it in items:
        s = str(it).strip()
        if len(s) == 1 and s.upper() in _SUPPLIER_CANON:
            codes.add(s.upper())
        else:
            codes |= _name_to_codes(s)
    return codes

# ====== primer3 (опционально) ======
_PRIMER3_AVAILABLE = False
_PRIMER3_USE_NEW = False
try:
    import primer3
    _PRIMER3_AVAILABLE = True
    _PRIMER3_USE_NEW = hasattr(primer3.bindings, "design_primers")
except Exception:
    _PRIMER3_AVAILABLE = False
    _PRIMER3_USE_NEW = False

# ====== I/O ======
Variant = namedtuple("Variant", "chrom pos ref alt mapped_id")

def load_names2(path: str) -> Dict[str, str]:
    """names2.txt: две колонки (input_chr -> fasta_id). Разделители: таб/пробел/запятая."""
    m = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[\t, ]+", line)
            if len(parts) < 2:
                continue
            key, val = parts[0], parts[1]
            m[key] = val
    return m

def detect_delim(header: str) -> str:
    if "\t" in header: return "\t"
    if "," in header: return ","
    # по умолчанию таб — ожидаемый формат
    return "\t"

def parse_input(path: str) -> List[Tuple[str,int,str,str]]:
    """Поддерживает .txt/.tsv (tab), CRLF/UTF-8/Windows-1251 — безопасно: errors='replace'."""
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        first = f.readline()
        if not first:
            return rows
        delim = detect_delim(first)
        hdr = [h.strip().lower() for h in first.strip("\r\n").split(delim)]
        # допускаем варианты имён колонок
        def idx(*cands):
            for c in cands:
                if c in hdr:
                    return hdr.index(c)
            return None

        i_chrom = idx("chrom", "chr", "chr1", "chromosome")
        i_pos   = idx("pos", "position")
        i_ref   = idx("ref", "reference")
        i_alt   = idx("alt", "alternate", "alt_allele")

        # если первая строка явно не заголовок, попробуем считать её как данные
        def parse_line(line: str):
            parts = line.strip("\r\n").split(delim)
            if len(parts) < 4:
                return None
            try:
                ch = parts[i_chrom] if i_chrom is not None else parts[0]
                po = int(parts[i_pos] if i_pos is not None else parts[1])
                rf = parts[i_ref] if i_ref is not None else parts[2]
                al = parts[i_alt] if i_alt is not None else parts[3]
                return (ch.replace("chr", "").strip(), po, rf.strip().upper(), al.strip().upper())
            except Exception:
                return None

        data = []
        # если похоже на заголовок (есть слова chrom/pos/ref/alt)
        if any(x in hdr for x in ("chrom","chr","pos","ref","alt","position","reference","alternate")):
            # читаем оставшиеся строки как данные
            for line in f:
                if not line.strip(): continue
                if "должно быть:" in line.lower():  # из ваших примеров
                    continue
                rec = parse_line(line)
                if rec: data.append(rec)
        else:
            # первая строка сразу данные
            rec = parse_line(first)
            if rec: data.append(rec)
            for line in f:
                if not line.strip(): continue
                rec = parse_line(line)
                if rec: data.append(rec)

        rows.extend(data)
    return rows

# ====== утилиты по ферментам ======
def cut_positions(enzyme, seq: str) -> List[int]:
    """Позиции разрезов (0..len) для линейной ДНК (включая 0 и конец)."""
    try:
        sites = enzyme.search(Seq(seq))  # 1-based позиции распознавания
    except Exception:
        sites = []
    cuts = [0]
    cuts += [p-1 for p in sites]  # 0-based
    L = len(seq)
    if L not in cuts:
        cuts.append(L)
    cuts = sorted(set(max(0, min(L, c)) for c in cuts))
    return cuts

def frag_lengths(cuts: List[int]) -> List[int]:
    return [cuts[i+1] - cuts[i] for i in range(len(cuts)-1)]

def frags_to_str(frags: List[int]) -> str:
    return "+".join(str(x) for x in frags) if frags else ""

def any_diag_delta(fr1: List[int], fr2: List[int], delta: int) -> bool:
    """Есть ли хотя бы одно различие фрагментов >= delta bp."""
    s1 = sorted(fr1)
    s2 = sorted(fr2)
    if len(s1) != len(s2):
        return True
    return any(abs(a-b) >= delta for a,b in zip(s1, s2))

def enzyme_site_string(enzyme) -> str:
    site = getattr(enzyme, "site", None)
    if isinstance(site, str):
        return site
    try:
        return str(site) if site is not None else ""
    except Exception:
        return ""

def enzyme_suppliers_human(enzyme) -> str:
    names = set()
    try:
        for v in enzyme.supplier_list():
            if v: names.add(str(v).strip())
    except Exception:
        pass
    if not names:
        codes = enzyme_supplier_codes(enzyme)
        for c in sorted(codes):
            aliases = _SUPPLIER_CANON.get(c, [])
            if aliases:
                nm = aliases[0]
                names.add(" ".join(w.capitalize() for w in nm.split()))
            else:
                names.add(c)
    return ", ".join(sorted(names))

# ====== primer3 ======
def primer3_pick(amplicon: str, snp_off: int, prod_min: int, prod_max: int, tm_min: float, tm_max: float) -> Optional[dict]:
    if not _PRIMER3_AVAILABLE:
        return None
    seq_args = {"SEQUENCE_TEMPLATE": amplicon, "SEQUENCE_TARGET": [snp_off, 1]}
    global_args = {
        "PRIMER_TASK": "pick_pcr_primers",
        "PRIMER_NUM_RETURN": 1,
        "PRIMER_MIN_TM": tm_min,
        "PRIMER_MAX_TM": tm_max,
        "PRIMER_OPT_TM": (tm_min + tm_max)/2.0,
        "PRIMER_PRODUCT_SIZE_RANGE": [[prod_min, prod_max]],
        "PRIMER_EXPLAIN_FLAG": 1,
    }
    try:
        if _PRIMER3_USE_NEW:
            res = primer3.bindings.design_primers(seq_args, global_args)
        else:
            res = primer3.bindings.designPrimers(seq_args, global_args)
        if res.get("PRIMER_PAIR_NUM_RETURNED", 0) > 0:
            return res
        return None
    except Exception:
        return None

# ====== основной расчёт ======
def main():
    p = argparse.ArgumentParser(description="Pick RFLP candidates with supplier filter (txt/tsv; no auto-reverse).")
    p.add_argument("--fasta", required=True, help="FASTA (indexed .fai)")
    p.add_argument("--input", required=True, help="variants txt/tsv; columns: chrom/chr, pos, ref, alt")
    p.add_argument("--names2", required=True, help="mapping of input chrom -> FASTA contig (two columns)")
    p.add_argument("--out", default="rflp_candidates.csv")
    p.add_argument("--flank", type=int, default=300)
    p.add_argument("--min-frag", type=int, default=80)
    p.add_argument("--max-frag", type=int, default=800)
    p.add_argument("--delta", type=int, default=25, help="min diagnostic delta between fragment sizes")
    p.add_argument("--gain-loss-only", action="store_true", help="require different number of cuts (gain/loss)")
    p.add_argument("--max-cuts", type=int, default=3, help="max cuts allowed per allele within amplicon")
    p.add_argument("--use-primer3", action="store_true")
    p.add_argument("--prod-min", type=int, default=250)
    p.add_argument("--prod-max", type=int, default=600)
    p.add_argument("--tm-min", type=float, default=58.0)
    p.add_argument("--tm-max", type=float, default=62.0)
    p.add_argument("--suppliers", default="", help="e.g. 'NEI' or 'NEB,Thermo,SibEnzyme'")
    args = p.parse_args()

    variants_raw = parse_input(args.input)
    if not variants_raw:
        print("[ERROR] input variants empty or unreadable", file=sys.stderr)
        sys.exit(2)
    names2 = load_names2(args.names2)

    allowed_codes = parse_suppliers_arg(args.suppliers)
    print(f"[INFO] Variants parsed: {len(variants_raw)}")
    if allowed_codes:
        kept = 0
        for e in AllEnzymes:
            if enzyme_supplier_codes(e) & allowed_codes:
                kept += 1
        humans = []
        for c in sorted(allowed_codes):
            aliases = _SUPPLIER_CANON.get(c, [])
            humans.append(aliases[0].title() if aliases else c)
        print(f"[INFO] Supplier filter: {kept}/{len(AllEnzymes)} enzymes kept for {humans}")
    else:
        print(f"[INFO] Supplier filter: not set (all enzymes allowed)")

    fa = Fasta(args.fasta, as_raw=True)

    out_rows = []
    header = [
        "variant","mapped_id","enzyme","site","pattern",
        "frags_ref","frags_alt","diag_delta_bp",
        "amplicon_start","amplicon_end","amplicon_len","snp_offset_in_amplicon",
        "primers_left","primers_right","tm_left","tm_right",
        "product_size","primer3_size","primer_left_start","primer_left_len","primer_right_start","primer_right_len",
        "suppliers"
    ]

    for chrom, pos, ref, alt in variants_raw:
        mapped_id = names2.get(chrom, chrom)
        start = max(1, pos - args.flank)
        end   = pos + args.flank
        try:
            seq_win = str(fa[mapped_id][start-1:end])
        except Exception as e:
            print(f"[WARN] Cannot fetch {mapped_id}:{start}-{end}: {e}", file=sys.stderr)
            continue

        snp_off = pos - start
        if snp_off < 0 or snp_off >= len(seq_win):
            print(f"[WARN] SNP offset out of window for {chrom}:{pos}", file=sys.stderr)
            continue

        genome_base = seq_win[snp_off].upper()
        if genome_base not in (ref.upper(), alt.upper()):
            print(f"[WARN] Reference base mismatch at {chrom}:{pos} FASTA has {genome_base} but input REF is {ref}", file=sys.stderr)

        seq_ref_win = (seq_win[:snp_off] + ref + seq_win[snp_off+1:]).upper()
        seq_alt_win = (seq_win[:snp_off] + alt + seq_win[snp_off+1:]).upper()

        pres = None
        amp_start0 = 0
        amp_end0_excl = len(seq_ref_win)
        amplicon_len = len(seq_ref_win)
        snp_off_amp = snp_off
        primers_left = primers_right = ""
        tm_left: Optional[float] = None
        tm_right: Optional[float] = None
        p3_size: Optional[int] = None
        pl_start = pl_len = pr_start = pr_len = None

        if args.use_primer3:
            pres = primer3_pick(seq_ref_win, snp_off, args.prod_min, args.prod_max, args.tm_min, args.tm_max)
            if pres:
                primers_left  = pres["PRIMER_LEFT_0_SEQUENCE"]
                primers_right = pres["PRIMER_RIGHT_0_SEQUENCE"]
                # numeric Tm
                tm_left = float(pres["PRIMER_LEFT_0_TM"])
                tm_right = float(pres["PRIMER_RIGHT_0_TM"])
                pl_start, pl_len = pres["PRIMER_LEFT_0"]
                pr_start, pr_len = pres["PRIMER_RIGHT_0"]
                p3_size = int(pres["PRIMER_PAIR_0_PRODUCT_SIZE"])

                amp_start0 = int(pl_start)
                amp_end0_excl = amp_start0 + int(pres["PRIMER_PAIR_0_PRODUCT_SIZE"])
                amp_start0 = max(0, min(len(seq_ref_win), amp_start0))
                amp_end0_excl = max(amp_start0, min(len(seq_ref_win), amp_end0_excl))

                amplicon_len = amp_end0_excl - amp_start0
                snp_off_amp = snp_off - amp_start0

        seq_ref_amp = seq_ref_win[amp_start0:amp_end0_excl]
        seq_alt_amp = seq_alt_win[amp_start0:amp_end0_excl]

        amplicon_start_genome = start + amp_start0
        amplicon_end_genome   = start + amp_end0_excl - 1

        enzymes = []
        for e in AllEnzymes:
            if allowed_codes and not (enzyme_supplier_codes(e) & allowed_codes):
                continue
            enzymes.append(e)

        hits_here = 0
        for enz in enzymes:
            cuts_r = cut_positions(enz, seq_ref_amp)
            cuts_a = cut_positions(enz, seq_alt_amp)

            if len(cuts_r)-2 > args.max_cuts or len(cuts_a)-2 > args.max_cuts:
                continue

            fr_r = frag_lengths(cuts_r)
            fr_a = frag_lengths(cuts_a)

            if (fr_r and (min(fr_r) < args.min_frag or max(fr_r) > args.max_frag)) or \
               (fr_a and (min(fr_a) < args.min_frag or max(fr_a) > args.max_frag)):
                continue

            if args.gain_loss_only and (len(cuts_r) == len(cuts_a)):
                continue

            if not any_diag_delta(fr_r, fr_a, args.delta):
                continue

            patt = "gain/loss" if (len(cuts_r) != len(cuts_a)) else "shift"
            site = enzyme_site_string(enz)

            out_rows.append([
                f"{chrom}:{pos} {ref}>{alt}",
                mapped_id,
                enz.__name__,
                site,
                patt,
                frags_to_str(fr_r),               # text
                frags_to_str(fr_a),               # text
                int(args.delta),                  # number
                int(amplicon_start_genome),       # number
                int(amplicon_end_genome),         # number
                int(amplicon_len),                # number
                int(snp_off_amp),                 # number
                primers_left,                     # text
                primers_right,                    # text
                (round(tm_left, 3) if tm_left is not None else None),   # number
                (round(tm_right, 3) if tm_right is not None else None), # number
                int(amplicon_len),                # product_size as number
                (int(p3_size) if p3_size is not None else None),        # number or blank
                (int(pl_start) if pl_start is not None else None),
                (int(pl_len) if pl_len is not None else None),
                (int(pr_start) if pr_start is not None else None),
                (int(pr_len) if pr_len is not None else None),
                enzyme_suppliers_human(enz),
            ])
            hits_here += 1

        if hits_here == 0:
            out_rows.append([
                f"{chrom}:{pos} {ref}>{alt}",
                mapped_id,
                "",
                "",
                "no_enzyme_found",
                "",
                "",
                int(args.delta),
                int(amplicon_start_genome),
                int(amplicon_end_genome),
                int(amplicon_len),
                int(snp_off_amp),
                primers_left,
                primers_right,
                (round(tm_left, 3) if tm_left is not None else None),
                (round(tm_right, 3) if tm_right is not None else None),
                int(amplicon_len),
                (int(p3_size) if p3_size is not None else None),
                (int(pl_start) if pl_start is not None else None),
                (int(pl_len) if pl_len is not None else None),
                (int(pr_start) if pr_start is not None else None),
                (int(pr_len) if pr_len is not None else None),
                "",
            ])

    out_lower = args.out.lower()
    if out_lower.endswith(".xlsx"):
        try:
            from openpyxl import Workbook
        except Exception:
            print("[ERROR] openpyxl is required to write .xlsx. Install with: pip install openpyxl", file=sys.stderr)
            sys.exit(3)
        wb = Workbook()
        ws = wb.active
        ws.title = "rflp_candidates"
        ws.append(header)
        for row in out_rows:
            ws.append(row)
        # Auto-width (basic): broaden numeric columns a bit
        for col_idx, h in enumerate(header, start=1):
            col_letter = ws.cell(row=1, column=col_idx).column_letter
            ws.column_dimensions[col_letter].width = max(12, min(60, len(str(h)) + 2))
        wb.save(args.out)
        print(f"Done: {args.out} ({len(out_rows)} rows)")
    else:
        with open(args.out, "w", newline="", encoding="utf-8") as w:
            writer = csv.writer(w)
            writer.writerow(header)
            writer.writerows(out_rows)
        print(f"Done: {args.out} ({len(out_rows)} rows)")

if __name__ == "__main__":
    main()
