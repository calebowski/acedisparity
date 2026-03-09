#!/bin/bash

JOB_FILE="job_ids.txt"

if [ ! -f "$JOB_FILE" ]; then
    echo "Error: Could not find $JOB_FILE."
    exit 1
fi

echo "Querying SLURM database for CPU Utilized..."
echo "------------------------------------------------------------------"
printf "%-12s | %-15s | %s\n" "Job ID" "Tasks Counted" "Total CPU (Hours)"
echo "------------------------------------------------------------------"

while IFS= read -r JOB_ID || [ -n "$JOB_ID" ]; do
    # Skip empty lines and strip hidden Windows carriage returns
    [ -z "$JOB_ID" ] && continue
    JOB_ID=$(echo "$JOB_ID" | tr -d '\r')

    # Query SLURM for JobID and Time
    RAW_DATA=$(sacct -j "$JOB_ID" -o JobID,TotalCPU -n 2>/dev/null)

    if [ -z "$RAW_DATA" ]; then
        printf "%-12s | %-15s | %s\n" "$JOB_ID" "0" "Not found"
        continue
    fi

    # Parse the output
    echo "$RAW_DATA" | awk -v jid="$JOB_ID" '
    BEGIN { array_secs=0; base_secs=0; array_count=0; base_count=0 }
    {
        # Skip the duplicate .b+, .batch, or .extern lines
        if ($1 ~ /\./) next;
        
        # Skip blank lines
        if ($2 == "") next;

        time_str = $2;
        
        # Strip off the milliseconds if they exist
        sub(/\..*$/, "", time_str);
        
        days=0;
        # Check if SLURM added days (DD-HH:MM:SS)
        if (time_str ~ /-/) {
            split(time_str, d, "-");
            days = d[1];
            time_str = d[2];
        }
        
        # Split the remaining time by colons
        n = split(time_str, t, ":");
        line_secs = 0;
        
        if (n == 3) {
            line_secs = days*86400 + t[1]*3600 + t[2]*60 + t[3];
        } else if (n == 2) {
            line_secs = days*86400 + t[1]*60 + t[2];
        }

        # SORTING LOGIC: Is this the Master row or an Array Task?
        if ($1 == jid) {
            # This is the Master summary row (or a single non-array job)
            base_secs = line_secs;
            base_count = 1;
        } else if ($1 ~ jid "_") {
            # This is an individual array task (e.g., 8558401_1)
            array_secs += line_secs;
            array_count++;
        }
    }
    END {
        # If we found array tasks, print their sum and ignore the Master row
        if (array_count > 0) {
            printf "%-12s | %-15d | %.2f\n", jid, array_count, array_secs/3600;
        } 
        # If we only found a Master row, print that (it was a single job)
        else if (base_count > 0) {
            printf "%-12s | %-15d | %.2f\n", jid, base_count, base_secs/3600;
        } 
        # If nothing valid was found
        else {
            printf "%-12s | %-15d | %s\n", jid, 0, "No valid data";
        }
    }'
done < "$JOB_FILE"

echo "------------------------------------------------------------------"