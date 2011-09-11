#!/bin/bash
for file in *.f *.F ; do
    echo UNVersioning file $file
    sed -e "/c___c/,/c___c/ d" $file > temp.tmp
    mv temp.tmp $file
done
