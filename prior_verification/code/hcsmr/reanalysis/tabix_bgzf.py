"""Minimal read-only Tabix/BGZF region reader for locked local files.

It exists solely as a portable fallback when a persistent htslib executable is
unavailable.  It supports standard ``TBI\x01`` indexes and reads bytes from a
local, already verified BGZF file.  It does not normalize or alter records.
"""
from __future__ import annotations

import gzip
import struct
from pathlib import Path


def reg2bins(start: int, end: int) -> list[int]:
    """Return Tabix bins overlapping a zero-based half-open interval."""
    if start < 0 or end <= start:
        raise ValueError("region must be a non-empty zero-based half-open interval")
    end -= 1
    bins = [0]
    for level in range(1, 6):
        shift = 14 + 3 * (5 - level)
        offset = ((1 << (3 * level)) - 1) // 7
        bins.extend(range((start >> shift) + offset, (end >> shift) + offset + 1))
    return bins


def _read_exact(handle, count: int) -> bytes:
    data = handle.read(count)
    if len(data) != count:
        raise ValueError("truncated Tabix/BGZF input")
    return data


def read_tbi(path: str | Path) -> tuple[list[str], dict[str, dict[int, list[tuple[int, int]]]]]:
    """Read reference names and bin chunks from a standard Tabix index."""
    # The eQTL Catalogue distributes this index as a gzip-wrapped TBI stream;
    # standard uncompressed TBI files are also accepted by ``gzip.open`` only
    # when detected first, so choose the reader from its magic bytes.
    raw = Path(path).open("rb")
    magic = raw.read(2)
    raw.seek(0)
    handle = gzip.GzipFile(fileobj=raw, mode="rb") if magic == b"\x1f\x8b" else raw
    try:
        if _read_exact(handle, 4) != b"TBI\x01":
            raise ValueError("not a standard TBI\\x01 index")
        n_ref, _format, _col_seq, _col_beg, _col_end, _meta, _skip, names_length = struct.unpack(
            "<8i", _read_exact(handle, 32)
        )
        names = _read_exact(handle, names_length).rstrip(b"\0").decode("utf-8").split("\0")
        if len(names) != n_ref:
            raise ValueError("Tabix reference-name count does not match header")
        all_bins: dict[str, dict[int, list[tuple[int, int]]]] = {}
        for name in names:
            bins: dict[int, list[tuple[int, int]]] = {}
            for _ in range(struct.unpack("<I", _read_exact(handle, 4))[0]):
                bin_id, chunk_count = struct.unpack("<II", _read_exact(handle, 8))
                bins[bin_id] = [struct.unpack("<QQ", _read_exact(handle, 16)) for _ in range(chunk_count)]
            interval_count = struct.unpack("<I", _read_exact(handle, 4))[0]
            _read_exact(handle, 8 * interval_count)
            all_bins[name] = bins
    finally:
        handle.close()
        if not raw.closed:
            raw.close()
    return names, all_bins


def _merge_chunks(chunks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(chunks):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def _bgzf_block(handle, compressed_offset: int) -> tuple[bytes, int]:
    handle.seek(compressed_offset)
    fixed = _read_exact(handle, 12)
    if fixed[:4] != b"\x1f\x8b\x08\x04":
        raise ValueError(f"missing BGZF header at compressed offset {compressed_offset}")
    extra_length = struct.unpack("<H", fixed[10:12])[0]
    extra = _read_exact(handle, extra_length)
    cursor = 0
    block_size: int | None = None
    while cursor < len(extra):
        subfield_id = extra[cursor:cursor + 2]
        subfield_length = struct.unpack("<H", extra[cursor + 2:cursor + 4])[0]
        subfield = extra[cursor + 4:cursor + 4 + subfield_length]
        if subfield_id == b"BC" and subfield_length == 2:
            block_size = struct.unpack("<H", subfield)[0] + 1
            break
        cursor += 4 + subfield_length
    if block_size is None:
        raise ValueError("BGZF BC subfield not found")
    handle.seek(compressed_offset)
    return gzip.decompress(_read_exact(handle, block_size)), block_size


def fetch_region(path: str | Path, index_path: str | Path, chromosome: str, start_1based: int, end_1based: int) -> list[str]:
    """Return complete text records in the requested one-based inclusive region."""
    if start_1based <= 0 or end_1based < start_1based:
        raise ValueError("region must be one-based inclusive with a positive start")
    _, indexes = read_tbi(index_path)
    ref = str(chromosome).removeprefix("chr")
    available = indexes.get(ref) or indexes.get("chr" + ref)
    if available is None:
        raise ValueError(f"chromosome {chromosome} is absent from the Tabix index")
    chunks = [chunk for bin_id in reg2bins(start_1based - 1, end_1based) for chunk in available.get(bin_id, [])]
    lines: list[str] = []
    with Path(path).open("rb") as handle:
        for chunk_start, chunk_end in _merge_chunks(chunks):
            compressed_offset, uncompressed_offset = chunk_start >> 16, chunk_start & 0xFFFF
            fragment = b""
            first_block = True
            while (compressed_offset << 16) < chunk_end:
                decoded, block_size = _bgzf_block(handle, compressed_offset)
                begin = uncompressed_offset if first_block else 0
                finish = len(decoded)
                if (compressed_offset << 16) == (chunk_end & ~0xFFFF):
                    finish = min(finish, chunk_end & 0xFFFF)
                fragment += decoded[begin:finish]
                if (compressed_offset << 16) == (chunk_end & ~0xFFFF):
                    break
                compressed_offset += block_size
                first_block = False
            lines.extend(line.decode("utf-8") for line in fragment.splitlines() if line)
    return lines
