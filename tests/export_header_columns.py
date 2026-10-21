"""Export schema names for native R tests without an additional R dependency."""

import csv
import json
import sys


def export_columns(source, destination):
    with open(source) as stream:
        files = json.load(stream)["files"]
    with open(destination, "w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["kind", "variant", "column_index", "column_name"])
        for kind, entry in files.items():
            for variant, schema in enumerate(entry["schemas"], 1):
                for index, name in enumerate(schema["columns"], 1):
                    writer.writerow([kind, variant, index, name])


if __name__ == "__main__":
    export_columns(*sys.argv[1:])
