#!/bin/bash

for file in *.fq.gz; do
    if ! gzip -t "$file"; then
        echo "❌ $file is corrupted or improperly compressed."
    else
        echo "✅ $file is valid."
    fi
done
