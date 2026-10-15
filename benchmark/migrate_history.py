"""Migrate dashboard data.js, downloaded JSON, or BenchmarkTools JSON without dropping samples."""

import argparse
import copy
import json
from pathlib import Path


RENAMES = json.loads(Path(__file__).with_name("hierarchy.json").read_text())
PREFIX = "window.BENCHMARK_DATA = "


def canonical_name(name):
    parts = RENAMES.get(name, name.split("/"))
    if len(parts) != 3 or any(not p or "/" in p for p in parts):
        raise ValueError(f"Unmapped benchmark name: {name!r}")
    return "/".join(parts)


def migrate_dashboard(data):
    for runs in data["entries"].values():
        for run in runs:
            names = set()
            for bench in run["benches"]:
                name = canonical_name(bench["name"])
                if name in names:
                    raise ValueError(f"Duplicate benchmark after migration: {name}")
                names.add(name)
                bench["name"] = name


def is_group(node):
    return isinstance(node, list) and len(node) == 2 and node[0] in (
        "BenchmarkGroup", "BenchmarkTools.BenchmarkGroup"
    )


def migrate_group(group, path=()):
    """Move leaves in place so existing group metadata and result objects survive."""
    data = group[1]["data"]
    for key, value in list(data.items()):
        child = (*path, key)
        if is_group(value):
            migrate_group(value, child)
            continue
        target = tuple(canonical_name("/".join(child)).split("/"))
        if child == target:
            continue
        # Historical drift only omitted the middle group under these topics.
        if len(path) != 1 or target[0] != path[0]:
            raise ValueError(f"Unsupported BenchmarkTools path: {child}")
        case, legend = target[1:]
        if case not in data:
            data[case] = [group[0], {"tags": [], "data": {}}]
        if not is_group(data[case]) or legend in data[case][1]["data"]:
            raise ValueError(f"Duplicate benchmark after migration: {'/'.join(target)}")
        data[case][1]["data"][legend] = value
        del data[key]


def migrate(data):
    result = copy.deepcopy(data)
    if isinstance(result, dict) and "entries" in result:
        migrate_dashboard(result)
    elif (isinstance(result, list) and len(result) == 2
          and isinstance(result[0], dict) and isinstance(result[1], list)
          and all(is_group(group) for group in result[1])):
        for group in result[1]:
            migrate_group(group)
    else:
        raise ValueError("Expected dashboard history or BenchmarkTools.save JSON")
    return result


def migrate_file(source, destination):
    text = source.read_text()
    javascript = text.startswith(PREFIX)
    payload = text[len(PREFIX):].strip().removesuffix(";") if javascript else text
    original = json.loads(payload)
    migrated = migrate(original)
    # Validate a second pass before touching the destination.
    if migrate(migrated) != migrated:
        raise ValueError("Migration is not idempotent")
    if migrated == original:
        if source != destination:
            destination.write_text(text)
        return False
    if source.resolve() == destination.resolve():
        backup = source.with_name(source.name + ".bak")
        if backup.exists():
            raise FileExistsError(f"Refusing to overwrite backup: {backup}")
        backup.write_bytes(source.read_bytes())
    rendered = json.dumps(migrated, ensure_ascii=False, indent=2) + "\n"
    destination.write_text(PREFIX + rendered if javascript else rendered)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path, nargs="?",
                        help="Omit for in-place migration with a .bak backup")
    args = parser.parse_args()
    changed = migrate_file(args.source, args.destination or args.source)
    print("Migrated" if changed else "Already canonical")


if __name__ == "__main__":
    main()
