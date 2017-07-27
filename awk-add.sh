#!/bin/awk
{ sum += $2 }
END { if (NR > 0) print sum / NR }
