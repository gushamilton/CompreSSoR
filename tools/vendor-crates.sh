#!/bin/sh
# Regenerate src/pcodec_native/vendor.tar.xz from src/pcodec_native/Cargo.lock.
# Needs network (and cargo >= the rust-version in Cargo.toml). Run after any
# change to Cargo.toml / Cargo.lock and commit the resulting tarball.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
crate="$root/src/pcodec_native"
work=$(mktemp -d "${TMPDIR:-/tmp}/compressor-vendor.XXXXXX")
trap 'rm -rf -- "$work"' EXIT HUP INT TERM

cargo vendor --locked --versioned-dirs --manifest-path "$crate/Cargo.toml" \
  "$work/vendor" >/dev/null
# Deterministic archive: sorted names, fixed owner and mtime.
(cd "$work" && find vendor -print | LC_ALL=C sort > list.txt)
if tar --version 2>/dev/null | grep -qi 'gnu tar'; then
  (cd "$work" && tar --no-recursion --owner=0 --group=0 --numeric-owner \
     --mtime='2000-01-01 00:00:00Z' -T list.txt -cf vendor.tar)
else
  (cd "$work" && tar -n --uid 0 --gid 0 --numeric-owner -T list.txt -cf vendor.tar)
fi
xz -9e -c "$work/vendor.tar" > "$crate/vendor.tar.xz"
ls -l "$crate/vendor.tar.xz"
