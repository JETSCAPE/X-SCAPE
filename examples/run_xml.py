#!/usr/bin/env python3
"""
examples/run_xml.py -- the two XMLs of one run_in_workdir.sh run, ready for its run directory.

    run_xml.py -d DATA_DIR -m MAIN_XML -u USER_XML -w RUN_DIR [--no-check]

Writes RUN_DIR/jetscape_main.xml and RUN_DIR/<user XML stem>.run.xml (the user's file is
not changed) and prints their paths, one per line. In both copies every element text that
starts with '../' is made absolute. Such paths are relative to the build directory (DATA_DIR),
the working directory of an in-build run, and not to RUN_DIR.

A '../' path with no file behind it is looked up in DATA_DIR by its trailing parts. An
installed tree (cmake --install, the macOS wheel's share/xscape) keeps iSS_tables,
iSS_parameters.dat and config/ directly in DATA_DIR, so
'../external_packages/iSS/iSS_tables' resolves to DATA_DIR/iSS_tables and
'../config/HH_LightExcitedMesons.xml' to DATA_DIR/config/HH_LightExcitedMesons.xml. The
names must match exactly, case included (macOS file systems ignore case).

Checks, before anything runs (exit status 3 when one fails, unless --no-check):

* every tag of the user XML is in the main XML, as runJetscape requires. It would stop at
  the first unknown tag with exit status 255; here every unknown tag is listed.
* every '../' path that the run uses has a file behind it. A path is used when the user XML
  has the module block it belongs to: <Stage><Module>, or the element itself for a
  <Stage><element> such as <IS><initial_profile_path>.

Only the standard library: X-SCAPE's build already needs Python 3.
"""

from __future__ import annotations

import argparse
import os
import sys
import xml.etree.ElementTree as ET

EXIT_CHECK = 3


def _exists_exact(path: str) -> bool:
    """``path`` exists with exactly this spelling in its last part. A case-insensitive file
    system would otherwise match MUSIC_input to music_input."""
    if not os.path.exists(path):
        return False
    parent, name = os.path.split(path.rstrip("/"))
    try:
        return name in os.listdir(parent or ".")
    except OSError:
        return False


def resolve(text: str, data_dir: str) -> tuple[str, bool]:
    """A '../' path -> (absolute path, whether it exists).

    First the path itself, relative to ``data_dir``. If that's missing, ``data_dir`` joined
    with each tail of its parts, longest first."""
    path = os.path.normpath(os.path.join(data_dir, text))
    if _exists_exact(path):
        return path, True
    parts = [p for p in os.path.normpath(text).split(os.sep) if p not in ("", "..", ".")]
    for i in range(len(parts)):
        cand = os.path.join(data_dir, *parts[i:])
        if _exists_exact(cand):
            return cand, True
    return path, False


def _elements(root):
    """(element, its tag path below the root) for every element below ``root``."""
    stack = [(child, (child.tag,)) for child in reversed(list(root))]
    while stack:
        el, path = stack.pop()
        yield el, path
        stack.extend((c, path + (c.tag,)) for c in reversed(list(el)))


def _tag_paths(root) -> set:
    return {path for el, path in _elements(root) if isinstance(el.tag, str)}


def _used(path: tuple, user_paths: set) -> bool:
    """A main-XML element at ``path`` (<Stage><Module><...>) is read by the run: the user
    XML has its module block and doesn't set the element itself (X-SCAPE reads the user's
    value first). Elements directly below a stage are not checked: either they select a
    module themselves (<IS><initial_profile_path>), and then the user XML sets them, or
    they are the stage's own settings."""
    return len(path) >= 3 and path[:2] in user_paths and path not in user_paths


def absolutize(tree, data_dir: str, user_paths: set | None) -> list[str]:
    """Make the tree's '../' texts absolute -> "<tag path> = <text>" of the missing ones.

    ``user_paths``: for the main XML, the user XML's tag paths, so that only the blocks the
    run uses are checked. None checks every path (the user XML itself)."""
    missing = []
    for el, path in _elements(tree.getroot()):
        text = (el.text or "").strip()
        if not text.startswith("../"):
            continue
        absolute, ok = resolve(text, data_dir)
        el.text = absolute
        if not ok and (user_paths is None or _used(path, user_paths)):
            missing.append(f"<{'><'.join(path)}> {text}")
    return missing


def unknown_tags(user_root, main_root) -> list[str]:
    """Tags of the user XML that the main XML lacks, as JetScape::CompareElementsFromXML
    checks them: by name, anywhere in the file."""
    known = {el.tag for el in main_root.iter()}
    seen, out = set(), []
    for el, path in _elements(user_root):
        if isinstance(el.tag, str) and el.tag not in known and el.tag not in seen:
            seen.add(el.tag)
            out.append(f"<{el.tag}> (in <{'><'.join(path)}>)")
    return out


def prepare(data_dir: str, main_xml: str, user_xml: str, run_dir: str,
            check: bool = True) -> tuple[str, str, list[str]]:
    """Write the two copies -> (main copy, user copy, problems found)."""
    main = ET.parse(main_xml)
    user = ET.parse(user_xml)
    problems = []
    unknown = unknown_tags(user.getroot(), main.getroot())
    if unknown:
        problems.append(
            f"tags of {user_xml} that the main XML {main_xml} lacks (runJetscape stops at "
            "the first: 'User XML tag <...> is unrecognized'). Either the tag is misspelled, "
            "or the XML is newer than this X-SCAPE: use an X-SCAPE whose "
            "config/jetscape_main.xml has the tag, or remove it:\n    "
            + "\n    ".join(unknown))
    missing = absolutize(user, data_dir, None)
    missing += absolutize(main, data_dir, _tag_paths(user.getroot()))
    if missing:
        problems.append(
            f"inputs with no file behind them, neither relative to {data_dir} nor in it. "
            "This X-SCAPE does not have them (a module or test data not built or not "
            "installed), or the path is wrong. Give your own inputs as absolute paths:\n    "
            + "\n    ".join(missing))
    os.makedirs(run_dir, exist_ok=True)
    main_out = os.path.join(run_dir, "jetscape_main.xml")
    stem = os.path.splitext(os.path.basename(user_xml))[0]
    user_out = os.path.join(run_dir, f"{stem}.run.xml")
    if os.path.abspath(user_out) == os.path.abspath(user_xml):
        user_out = os.path.join(run_dir, f"{stem}.run.run.xml")
    main.write(main_out, encoding="utf-8", xml_declaration=True)
    user.write(user_out, encoding="utf-8", xml_declaration=True)
    return main_out, user_out, problems if check else []


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="run_xml.py",
        description="The XMLs of one run_in_workdir.sh run: '../' paths made absolute, "
                    "the user XML's tags and the run's inputs checked.")
    p.add_argument("-d", dest="data_dir", required=True, help="the build or installed tree")
    p.add_argument("-m", dest="main_xml", required=True)
    p.add_argument("-u", dest="user_xml", required=True)
    p.add_argument("-w", dest="run_dir", required=True)
    p.add_argument("--no-check", action="store_true",
                   help="write the copies even if a check fails, and exit 0")
    a = p.parse_args(argv)
    try:
        main_out, user_out, problems = prepare(os.path.abspath(a.data_dir), a.main_xml,
                                               a.user_xml, a.run_dir, check=not a.no_check)
    except (OSError, ET.ParseError) as e:
        print(f"run_xml.py: {e}", file=sys.stderr)
        return EXIT_CHECK
    if problems:
        for msg in problems:
            print(f"ERROR: {msg}", file=sys.stderr)
        print("Nothing was run. (run_in_workdir.sh -f runs anyway.)", file=sys.stderr)
        return EXIT_CHECK
    print(main_out)
    print(user_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
