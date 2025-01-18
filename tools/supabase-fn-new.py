from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

from dotenv import find_dotenv, load_dotenv
from supabase import create_client

if len(sys.argv) < 2:
    sys.exit("Error: missing function name argument")

load_dotenv(find_dotenv())


if not os.environ.get("SUPABASE_SR_URL"):
    sys.exit("Error: missing SUPABASE_SR_URL environment variable")
if not os.environ.get("SUPABASE_SR_SERVICE_ROLE_KEY"):
    sys.exit("Error: missing SUPABASE_SR_SERVICE_ROLE_KEY environment variable")

supabase = create_client(
    supabase_url=os.environ["SUPABASE_SR_URL"],
    supabase_key=os.environ["SUPABASE_SR_SERVICE_ROLE_KEY"],
)

fn_name = sys.argv[1]
fn_path = Path(f"supabase/functions/{fn_name}")

subprocess.run(["supabase", "functions", "new", fn_name], check=False)  # noqa: S603, S607

import_map_json = fn_path / "import_map.json"
packages = sys.argv[2:]
import_dict = {"imports": {}}
pkg_dict = {"{pkg}": "https://esm.sh/{pkg}@{version}/"}
# @langchain/core@0.0.0/foo/bar -> @langchain/core
pkg_re = re.compile(r"^@([^/@]+)((/[^/@]+)+)?(/)?")


def get_npm_package_version(pkg: str):  # noqa: ANN201
    if not pkg.endswith("/"):
        print(f"WARNING: / suffix required in {pkg} if you need to access submodules")
    try:
        pkg = pkg_re.match(pkg).groups()
        result = subprocess.run(  # noqa: S603
            ["npm", "show", pkg, "version"],  # noqa: S607
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.decode("utf-8").strip()
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        return None


if packages:
    for pkg in packages:
        # parse to exxpected format
        pkg_cleaned = pkg_re.match(pkg)
        if not pkg_cleaned:
            print(f"WARNING: cannot parse {pkg}, skipping ...")
            continue
        pkg_cleaned = pkg_cleaned.group()
        version = (
            subprocess.run(
                ["npm", "show", pkg_cleaned, "version"],
                capture_output=True,
                check=False,
            )
            .stdout.decode("utf-8")
            .strip()
        )
        if not re.match(r"^\d+\.\d+\.\d+$", version):
            print(
                f"WARNING: cannot parse version {version} for {pkg_cleaned}, skipping ..."
            )
            continue
        import_pkg_dict = {
            k.format(pkg=pkg_cleaned): v.format(pkg=pkg_cleaned, version=version)
            for k, v in pkg_dict.items()
        }
        import_dict["imports"][pkg_cleaned] = import_pkg_dict[pkg_cleaned]
        print(
            f"Added package {pkg_cleaned} with version {version} to function {fn_name}"
        )

if import_dict["imports"]:
    import_map_json.write_text(json.dumps(import_dict, indent=2))
    print(f"Wrote {import_map_json}")
