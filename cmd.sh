#!/bin/bash
bedtools intersect -wao -a golga.merged.bed  -b $(cat asmbeds.txt) > overlap.by.asm.bed

