# queue

Utilities for interacting with [Reed-CompBio/spras-benchmarking-queue](https://github.com/Reed-CompBio/spras-benchmarking-queue).

This will be explained inside `/htcondor`, but the HTCondor cron job, every time it runs, internally checks against its lists of timestamps that it
has already ran.

The queue has files formatted as `{YYYY-MM}/{unix-timestamp}.txt`, where the `.txt` file contains
the `spras-benchmarking` commit at its time of creation.

`amend.py` auto-commits to `spras-benchmarking-queue` by cloning it locally.
