#!/usr/bin/env python3
"""Hold a per-profile advisory lock until the parent closes standard input."""
import argparse
import fcntl
import json
import os
import socket
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--path", required=True)
args = parser.parse_args()
flags = os.O_CREAT | os.O_RDWR
if hasattr(os, "O_NOFOLLOW"):
    flags |= os.O_NOFOLLOW
try:
    fd = os.open(args.path, flags, 0o600)
    fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
except (IOError, OSError):
    print("BUSY", flush=True)
    sys.exit(2)
try:
    os.fchmod(fd, 0o600)
    os.ftruncate(fd, 0)
    os.write(fd, json.dumps({"pid": os.getpid(), "host": socket.gethostname()}).encode("utf-8"))
    os.fsync(fd)
    print("LOCKED", flush=True)
    # EOF also releases the lock after a killed R parent; the lock filename
    # remaining on disk does not mean it is still locked.
    sys.stdin.buffer.read()
finally:
    os.close(fd)
