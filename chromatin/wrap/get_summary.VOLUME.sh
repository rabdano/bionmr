#!/bin/bash
grep -m 1000 VOLUME *.out | cut -c 4-8,78-92 > summary.VOLUME