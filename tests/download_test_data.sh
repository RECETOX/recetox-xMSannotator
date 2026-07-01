#!/bin/bash
#
# Download test data files from remote URLs defined in tests/remote-files/
# Each fetch_<folder>.txt file contains a list of URLs to download
#
# Usage:
#   ./download_test_data.sh              # Download all files from all fetch_*.txt files
#   ./download_test_data <folder>        # Download files for a specific folder (e.g., sourceforge)
#   ./download_test_data --dry-run       # Show what would be downloaded without downloading
#   ./download_test_data --list          # List all available folders and their files
#
# Files are downloaded to: tests/test-data/<folder>/
#

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOTE_FILES_DIR="${SCRIPT_DIR}/remote-files"
TEST_DATA_DIR="${SCRIPT_DIR}/testthat/test-data"

DRY_RUN=false
SPECIFIC_FOLDER=""

# Function to print status messages
log_info() {
    echo "[INFO] $1"
}

log_success() {
    echo "[SUCCESS] $1"
}

log_warning() {
    echo "[WARNING] $1"
}

log_error() {
    echo "[ERROR] $1"
}

# Show usage
show_usage() {
    echo "Usage: $0 [OPTIONS] [FOLDER]"
    echo ""
    echo "Options:"
    echo "  --dry-run    Show what would be downloaded without downloading"
    echo "  --list       List all available folders and their files"
    echo "  --help       Show this help message"
    echo ""
    echo "Arguments:"
    echo "  FOLDER       Specific folder to download (e.g., sourceforge, qc_solvent)"
    echo ""
    echo "Examples:"
    echo "  $0                      # Download all files from all folders"
    echo "  $0 sourceforge          # Download only sourceforge files"
    echo "  $0 --dry-run            # Preview downloads"
    echo "  $0 --list               # List available folders"
}

# List all available folders
list_folders() {
    log_info "Available folders:"
    echo ""
    for fetch_file in "${REMOTE_FILES_DIR}"/fetch_*.txt; do
        if [[ -f "$fetch_file" ]]; then
            folder_name=$(basename "$fetch_file" | sed 's/fetch_//' | sed 's/\.txt//')
            file_count=$(grep -c "^http" "$fetch_file" 2>/dev/null || echo "0")
            echo "  - ${folder_name} (${file_count} files)"
        fi
    done
    echo ""
}

# Get folder name from fetch file
get_folder_name() {
    local fetch_file="$1"
    basename "$fetch_file" | sed 's/fetch_//' | sed 's/\.txt//'
}

# Download files from a single fetch file
download_folder() {
    local fetch_file="$1"
    local folder_name=$(get_folder_name "$fetch_file")
    local target_dir="${TEST_DATA_DIR}/${folder_name}"

    log_info "Processing: ${folder_name}"

    # Create target directory
    if [[ ! -d "$target_dir" ]]; then
        if [[ "$DRY_RUN" == true ]]; then
            log_info "[DRY-RUN] Would create: ${target_dir}"
        else
            mkdir -p "$target_dir"
            log_info "Created directory: ${target_dir}"
        fi
    fi

    # Read URLs from fetch file and download
    local success_count=0
    local fail_count=0
    local line_num=0

    while IFS= read -r url || [[ -n "$url" ]]; do
        line_num=$((line_num + 1))

        # Skip empty lines and comments
        [[ -z "$url" ]] && continue
        [[ "$url" =~ ^[[:space:]]*# ]] && continue

        # Extract filename from URL
        local filename=$(basename "$url")
        local target_file="${target_dir}/${filename}"

        if [[ "$DRY_RUN" == true ]]; then
            log_info "[DRY-RUN] Would download: ${filename}"
            success_count=$((success_count + 1))
            continue
        fi

        log_info "  Downloading: ${filename}"

        # Download the file
        if curl -sL "$url" -o "${target_file}.tmp"; then
            # Check if download was successful (file has content)
            if [[ -s "${target_file}.tmp" ]]; then
                mv "${target_file}.tmp" "$target_file"
                local size=$(du -h "$target_file" | cut -f1)
                log_success "    Saved: ${filename} (${size})"
                success_count=$((success_count + 1))
            else
                rm -f "${target_file}.tmp"
                log_error "    Failed: ${filename} (empty download)"
                fail_count=$((fail_count + 1))
            fi
        else
            rm -f "${target_file}.tmp"
            log_error "    Failed: ${filename} (download error)"
            fail_count=$((fail_count + 1))
        fi
    done < "$fetch_file"

    log_info "  Completed: ${success_count} succeeded, ${fail_count} failed"
    echo ""

    return $fail_count
}

# Main execution
main() {
    echo "=========================================="
    echo "Test Data Download Script"
    echo "=========================================="
    echo ""

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            --list)
                list_folders
                exit 0
                ;;
            --help|-h)
                show_usage
                exit 0
                ;;
            *)
                SPECIFIC_FOLDER="$1"
                shift
                ;;
        esac
    done

    # Check if remote-files directory exists
    if [[ ! -d "$REMOTE_FILES_DIR" ]]; then
        log_error "Remote files directory not found: ${REMOTE_FILES_DIR}"
        exit 1
    fi

    # Show dry-run mode
    if [[ "$DRY_RUN" == true ]]; then
        log_info "Running in dry-run mode - no files will be downloaded"
        echo ""
    fi

    # Determine which folders to process
    local folders_to_process=()

    if [[ -n "$SPECIFIC_FOLDER" ]]; then
        # Process specific folder
        local fetch_file="${REMOTE_FILES_DIR}/fetch_${SPECIFIC_FOLDER}.txt"
        if [[ ! -f "$fetch_file" ]]; then
            log_error "Fetch file not found: ${fetch_file}"
            log_info "Available folders:"
            list_folders
            exit 1
        fi
        folders_to_process+=("$fetch_file")
    else
        # Process all folders
        for fetch_file in "${REMOTE_FILES_DIR}"/fetch_*.txt; do
            if [[ -f "$fetch_file" ]]; then
                folders_to_process+=("$fetch_file")
            fi
        done
    fi

    if [[ ${#folders_to_process[@]} -eq 0 ]]; then
        log_error "No fetch files found to process"
        exit 1
    fi

    # Process each folder
    local total_success=0
    local total_fail=0

    for fetch_file in "${folders_to_process[@]}"; do
        download_folder "$fetch_file"
        local result=$?
        total_fail=$((total_fail + result))
    done

    # Summary
    echo "=========================================="
    echo "Download Summary"
    echo "=========================================="
    if [[ "$DRY_RUN" == true ]]; then
        log_info "Dry-run complete - no files were downloaded"
    elif [[ $total_fail -eq 0 ]]; then
        log_success "All files downloaded successfully!"
    else
        log_error "Some downloads failed: ${total_fail}"
        exit 1
    fi
}

# Run main function
main "$@"
