#!/bin/bash

shopt -s expand_aliases

#
# A script to run predefined sequences of commands for specific numbers.
# It loads tasks from an external script file provided by the user.
#

# --- Usage function ---
# This function is displayed if the script is run with -h/--help or with incorrect arguments.
usage() {
    echo "Usage: $0 -t <tasks_file> [number_spec]"
    echo "Runs predefined command sequences for specific numbers by loading a task script."
    echo ""
    echo "Arguments:"
    echo "  [number_spec]   (Optional) A specification of which tasks to run. If omitted,"
    echo "                  all tasks defined in the task script will be executed. It can be:"
    echo "                  - a single number (e.g., 5)"
    echo "                  - a range (e.g., 1-10)"
    echo "                  - a comma-separated list (e.g., 1,3,5)"
    echo "                  - a mix of range and comma-separated list (e.g., 1-3,5)"
    echo ""
    echo "Options:"
    echo "  -t, --tasks     (Required) The script file containing task definitions."
    echo "  -h, --help      Display this help message."
    echo ""
    echo "Example:"
    echo "  # Run tasks 1-3 from the 'my_tasks.sh' file"
    echo "  $0 -t my_tasks 1-3"
    echo ""
    echo "  # Run all tasks from 'my_tasks.sh'"
    echo "  $0 -t my_tasks"
    exit 1
}

# --- Some globals ---
export SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export PATH=$SCRIPTDIR:$PATH
export EXPDIR=`dirname $SCRIPTDIR`
export CODEDIR=$EXPDIR/code
export LOGDIR=$EXPDIR/logs
export OUTDIR=$EXPDIR/out
export TMPDIR=$EXPDIR/tmp
export RESDIR=$EXPDIR/results

export PATH=$EXPDIR/code/minimap2/:$PATH
export PATH=$EXPDIR/code/collinearity/build/:$PATH
export PATH=$EXPDIR/code/rawhash2/:$PATH
export PATH=$EXPDIR/code/metagraph/metagraph/build/:$PATH
export SPUMONI_BUILD_DIR=$EXPDIR/code/spumoni/build/
export PATH=$SPUMONI_BUILD_DIR:$PATH
export TIMESTAMP=$(date +"%b%d_%H%M%S")
alias measure="/usr/bin/time -f \"CPU=%P\nElapsed=%E\nMaxRSS=%M KB\""

# --- Parse command-line arguments ---
TASKS_NAME=""
TASKS_FILE=""

# Loop to handle options first
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -t|--tasks)
        TASKS_NAME="$2"
        if [[ "$TASKS_NAME" == *.sh ]]; then
            TASKS_NAME="${TASKS_NAME%.sh}"
        fi
        TASKS_FILE="$SCRIPTDIR/$TASKS_NAME.sh"
        shift # past argument
        shift # past value
        ;;
        -h|--help)
        usage
        ;;
        *)    # This is not an option, so it must be the number specification
        if [ -z "$NUMBER_SPEC" ]; then
            NUMBER_SPEC="$1"
            shift # past argument
        else
            echo "Error: Unknown or duplicate argument '$1'" >&2
            usage
        fi
        ;;
    esac
done

# --- Validate that a tasks file was provided ---
if [ -z "$TASKS_FILE" ]; then
    echo "Error: No tasks file provided. The -t or --tasks option is required." >&2
    usage
fi

# --- Source the Task Script ---
if [ ! -f "$TASKS_FILE" ]; then
    echo "Error: Tasks file '$TASKS_FILE' not found." >&2
    exit 1
fi

export LOGFILE=$LOGDIR/${TASKS_NAME}_$TIMESTAMP.log

# --- Task Execution Function ---
# This function dynamically calls the appropriate task function based on the number.
run_task_for() {
    local number=$1
    local task_function="f$number"

    echo "--- Processing task: $number ---"

    # Check if a function with the constructed name exists.
    if declare -f "$task_function" > /dev/null; then
        # If it exists, call it.
        "$task_function"
    else
        # Otherwise, show a warning.
        echo "Warning: No task defined for number '$number'. Skipping."
    fi

    echo "" # Add a blank line for readability
}

# 'source' runs the script in the current shell, making its functions and variables available.
source "$TASKS_FILE"

# --- Main logic ---

# If no number specification is given, run all tasks defined in the sourced file.
if [ -z "$NUMBER_SPEC" ]; then
    # The ALL_TASKS variable is expected to be defined in the tasks file.
    if [ -z "$ALL_TASKS" ]; then
        echo "Warning: No tasks specified and the 'ALL_TASKS' variable is not defined in '$TASKS_NAME'."
        exit 0
    fi
    echo "No number specification provided. Running all defined tasks from '$TASKS_NAME'."
    NUMBER_SPEC="$ALL_TASKS"
fi

echo "Starting script with specification: '$NUMBER_SPEC' using tasks from '$TASKS_NAME'"
echo ""

# Use 'tr' to replace commas with spaces, allowing a 'for' loop to iterate over each part.
for item in $(echo "$NUMBER_SPEC" | tr ',' ' '); do
    # Check if the item is a range (e.g., "1-10").
    if [[ $item =~ ^[0-9]+-[0-9]+$ ]]; then
        start=$(echo "$item" | cut -d'-' -f1)
        end=$(echo "$item" | cut -d'-' -f2)

        if [ "$start" -gt "$end" ]; then
            echo "Warning: Invalid range '$item' (start > end). Skipping."
            continue
        fi
        
        for i in $(seq "$start" "$end"); do
            run_task_for "$i" 2>&1 | tee -a $LOGFILE
        done
    # Check if the item is a single number (e.g., "5").
    elif [[ $item =~ ^[0-9]+$ ]]; then
        run_task_for "$item" 2>&1 | tee -a $LOGFILE
    else
        echo "Warning: Invalid item '$item' in specification. Skipping."
    fi
done

echo "Script finished. Logs written to $LOGFILE"

