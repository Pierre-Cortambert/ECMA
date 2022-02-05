#!/bin/sh

algorithmes = "s"
instances = "20_USA-road-d.BAY.gr  20_USA-road-d.COL.gr  20_USA-road-d.NY.gr 
             40_USA-road-d.BAY.gr  40_USA-road-d.COL.gr  40_USA-road-d.NY.gr 
             60_USA-road-d.BAY.gr  60_USA-road-d.COL.gr  60_USA-road-d.NY.gr
             80_USA-road-d.BAY.gr  80_USA-road-d.COL.gr  80_USA-road-d.NY.gr
             100_USA-road-d.BAY.gr 100_USA-road-d.COL.gr 100_USA-road-d.NY.gr
             120_USA-road-d.BAY.gr 120_USA-road-d.COL.gr 120_USA-road-d.NY.gr 
             140_USA-road-d.BAY.gr 140_USA-road-d.COL.gr 140_USA-road-d.NY.gr 
             160_USA-road-d.BAY.gr 160_USA-road-d.COL.gr 160_USA-road-d.NY.gr
             180_USA-road-d.BAY.gr 180_USA-road-d.COL.gr 180_USA-road-d.NY.gr "

for algo in $algorithmes;do
    for inst in $instances; do
        julia main.jl --algo $algo --file $inst
    done
done

