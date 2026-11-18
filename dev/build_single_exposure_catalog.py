"""One-off: reduce exposures.csv to one row per single-exposure plate.

Drops plates with more than one distinct real exposure (exposureNumber >= 0 -
negative values look like "merged" summary rows over multiple exposures, not
separate exposures themselves), since a multi-exposure plate photo shows
overlapping star fields that a single WCS solve isn't well-defined for.

Run once: python dev/build_single_exposure_catalog.py
"""
import csv
from collections import defaultdict
from pathlib import Path

SRC = Path(__file__).parent / "exposures.csv"
DEST = Path(__file__).parent / "exposures_single.csv"

OUT_COLUMNS = [
    "plateid", "naxis1", "naxis2", "ctr_ra", "ctr_dec",
    "delta_ra_x", "delta_ra_y", "delta_dec_x", "delta_dec_y", "datetime",
]

DELETED_VALUES = {"1", "1.0"}


def main() -> None:
    rows_by_plate: dict[str, list[dict[str, str]]] = defaultdict(list)
    with open(SRC, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("deleted") in DELETED_VALUES:
                continue
            rows_by_plate[row["plateid"]].append(row)

    kept = 0
    dropped_multi = 0
    dropped_no_zero = 0

    with open(DEST, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=OUT_COLUMNS)
        writer.writeheader()
        for plate_id, rows in rows_by_plate.items():
            real_exposures = {r["exposureNumber"] for r in rows if int(r["exposureNumber"]) >= 0}
            if len(real_exposures) > 1:
                dropped_multi += 1
                continue
            zero_row = next((r for r in rows if r["exposureNumber"] == "0"), None)
            if zero_row is None:
                dropped_no_zero += 1
                continue
            writer.writerow({col: zero_row[col] for col in OUT_COLUMNS})
            kept += 1

    print(f"kept={kept} dropped_multi_exposure={dropped_multi} dropped_no_exposure_zero={dropped_no_zero}")
    print(f"wrote {DEST}")


if __name__ == "__main__":
    main()
