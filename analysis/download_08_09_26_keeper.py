#!/usr/bin/env python3
"""Download the user-supplied Keeper share and verify every JPK archive."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import zipfile

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


ROOT = Path(__file__).resolve().parents[1]
SHARE_TOKEN = "6f30d59535114f89b568"
HOST = "https://keeper.mpdl.mpg.de"
DESTINATION = ROOT / "raw" / f"keeper_{SHARE_TOKEN}"
DATA = DESTINATION / "08-09-26"


def session() -> requests.Session:
    client = requests.Session()
    retry = Retry(total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504])
    client.mount("https://", HTTPAdapter(max_retries=retry))
    return client


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inventory(client: requests.Session, folder: str = "/") -> list[dict]:
    response = client.get(
        f"{HOST}/api/v2.1/share-links/{SHARE_TOKEN}/dirents/",
        params={"path": folder}, timeout=(15, 60),
    )
    response.raise_for_status()
    payload = response.json()
    if folder == "/":
        (DESTINATION / "root_dirents.json").write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )
    files = []
    for item in payload["dirent_list"]:
        if item["is_dir"]:
            files.extend(inventory(client, item["folder_path"]))
        else:
            files.append(item)
    return files


def download(item: dict) -> dict:
    relative = PurePosixPath(item["file_path"].lstrip("/"))
    if not relative.parts or any(part in (".", "..") or ":" in part or "\\" in part for part in relative.parts):
        raise ValueError(f"Unsafe relative path: {relative}")
    path = DATA.joinpath(*relative.parts).resolve()
    path.relative_to(DATA.resolve())
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.stat().st_size != item["size"]:
        temporary = path.with_name(path.name + ".part")
        with session() as client, client.get(
            f"{HOST}/d/{SHARE_TOKEN}/files/",
            params={"p": item["file_path"], "dl": 1},
            stream=True, timeout=(15, 60),
        ) as response:
            response.raise_for_status()
            with temporary.open("wb") as stream:
                for chunk in response.iter_content(1024 * 1024):
                    if chunk:
                        stream.write(chunk)
        if temporary.stat().st_size != item["size"]:
            raise RuntimeError(f"Download size mismatch: {path.name}")
        temporary.replace(path)
    crc = None
    if path.suffix in (".jpk-force-map", ".jpk-force", ".jpk-qi-data"):
        with zipfile.ZipFile(path) as archive:
            crc = archive.testzip()
            if crc is not None:
                raise RuntimeError(f"Failed archive CRC: {path.name}, {crc}")
    return {
        "share_relative_path": item["file_path"],
        "local_path": path.relative_to(ROOT).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
        "remote_last_modified": item["last_modified"],
        "jpk_zip_crc": "passed" if crc is None and zipfile.is_zipfile(path) else "not_applicable",
    }


def main() -> None:
    DESTINATION.mkdir(parents=True, exist_ok=True)
    with session() as client:
        files = inventory(client)
    print(f"Share: {len(files)} files, {sum(item['size'] for item in files):,} bytes", flush=True)
    records = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(download, item) for item in files]
        for future in as_completed(futures):
            record = future.result()
            records.append(record)
            print(f"[{len(records):02d}/{len(files)}] verified {Path(record['local_path']).name}", flush=True)
    manifest = {
        "share_url": f"{HOST}/d/{SHARE_TOKEN}/",
        "download_verified_utc": datetime.now(timezone.utc).isoformat(),
        "total_files": len(records),
        "total_bytes": sum(record["bytes"] for record in records),
        "files": sorted(records, key=lambda item: item["local_path"]),
    }
    (DESTINATION / "download_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote {DESTINATION / 'download_manifest.json'}", flush=True)


if __name__ == "__main__":
    main()
