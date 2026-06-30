#!/bin/bash
#
# Download test data files from GitLab
# https://gitlab.ics.muni.cz/umsa/umsa-files/-/tree/master/testdata/recetox-xMSannotator
#
# Usage:
#   ./download_test_data.sh              # Download all files (requires GITLAB_TOKEN for private repo)
#   ./download_test_data --skip-auth     # Skip auth check and try public access
#   ./download_test_data --dry-run       # Show what would be downloaded without downloading
#
# Environment variables:
#   GITLAB_TOKEN - Personal access token for GitLab (required if repo is private)
#   GITLAB_PROJECT_ID - Project ID (default: auto-detected)
#

set -e

# Configuration
GITLAB_BASE="https://gitlab.ics.muni.cz"
PROJECT_PATH="umsa/umsa-files"
BRANCH="master"
BASE_DIR="testdata/recetox-xMSannotator"
TEST_DATA_DIR="$(dirname "$0")/testthat/test-data"

# Files to download organized by directory
declare -A FILES=(
    # sourceforge directory
    ["sourceforge/tempobjects.Rda"]="sourceforge/tempobjects.Rda"
    ["sourceforge/global_cor_integ.Rda"]="sourceforge/global_cor_integ.Rda"
    ["sourceforge/chemscoremat.Rds"]="sourceforge/chemscoremat.Rds"
    # qc_solvent directory
    ["qc_solvent/tempobjects.Rda"]="qc_solvent/tempobjects.Rda"
    ["qc_solvent/global_cor_integ.Rda"]="qc_solvent/global_cor_integ.Rda"
    ["qc_solvent/chemscoremat.Rds"]="qc_solvent/chemscoremat.Rds"
    # qc_matrix directory
    ["qc_matrix/tempobjects.Rda"]="qc_matrix/tempobjects.Rda"
    ["qc_matrix/global_cor_integ.Rda"]="qc_matrix/global_cor_integ.Rda"
    ["qc_matrix/chemscoremat.Rds"]="qc_matrix/chemscoremat.Rds"
    # batch1_neg directory
    ["batch1_neg/global_cor_integ.Rda"]="batch1_neg/global_cor_integ.Rda"
    ["batch1_neg/chemscoremat.Rds"]="batch1_neg/chemscoremat.Rds"
)

DRY_RUN=false
SKIP_AUTH_CHECK=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --skip-auth)
            SKIP_AUTH_CHECK=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--dry-run] [--skip-auth]"
            exit 1
            ;;
    esac
done

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

# Check if GitLab token is available
check_auth() {
    if [[ -n "$GITLAB_TOKEN" ]]; then
        log_info "GitLab token found, will use authenticated requests"
        return 0
    else
        log_warning "No GitLab token found. Set GITLAB_TOKEN environment variable for authenticated access."
        log_info "Example: export GITLAB_TOKEN=glpat-xxxxxxxxxxxxx"
        return 1
    fi
}

# Create directory structure
create_directories() {
    log_info "Creating directory structure..."
    for dir in sourceforge qc_solvent qc_matrix batch1_neg; do
        local full_path="$TEST_DATA_DIR/$dir"
        if [[ ! -d "$full_path" ]]; then
            if [[ "$DRY_RUN" == true ]]; then
                log_info "[DRY-RUN] Would create: $full_path"
            else
                mkdir -p "$full_path"
                log_info "Created: $full_path"
            fi
        else
            log_info "Directory exists: $full_path"
        fi
    done
}

# Download a single file
download_file() {
    local relative_path=$1
    local local_path="$TEST_DATA_DIR/$relative_path"
    local file_dir=$(dirname "$local_path")
    local file_name=$(basename "$relative_path")

    # Build the URL - try both raw and API endpoints
    local raw_url="$GITLAB_BASE/$PROJECT_PATH/-/raw/$BRANCH/$BASE_DIR/$relative_path"
    local api_url="$GITLAB_BASE/api/v4/projects/$PROJECT_PATH/repository/files/$BASE_DIR%2F$relative_path/raw?ref=$BRANCH"

    log_info "Downloading: $relative_path"

    if [[ "$DRY_RUN" == true ]]; then
        log_info "[DRY-RUN] Would download from: $raw_url"
        log_info "[DRY-RUN] Would save to: $local_path"
        return 0
    fi

    # Ensure directory exists
    mkdir -p "$file_dir"

    # Download with optional authentication
    if [[ -n "$GITLAB_TOKEN" ]]; then
        curl -sL \
            -H "PRIVATE-TOKEN: $GITLAB_TOKEN" \
            "$raw_url" \
            -o "$local_path.tmp"
    else
        curl -sL "$raw_url" -o "$local_path.tmp"
    fi

    # Check if download was successful
    if [[ -f "$local_path.tmp" ]] && [[ -s "$local_path.tmp" ]]; then
        mv "$local_path.tmp" "$local_path"
        local size=$(du -h "$local_path" | cut -f1)
        log_success "Downloaded: $file_name ($size)"
        return 0
    else
        rm -f "$local_path.tmp"
        log_error "Failed to download: $relative_path"
        return 1
    fi
}

# List all files that would be downloaded
list_files() {
    log_info "Files to download:"
    echo ""
    for key in "${!FILES[@]}"; do
        echo "  - $key"
    done | sort
    echo ""
    echo "Total: ${#FILES[@]} files"
}

# Main execution
main() {
    echo "=========================================="
    echo "Test Data Download Script"
    echo "=========================================="
    echo ""

    # Show what we're doing
    list_files

    if [[ "$DRY_RUN" == true ]]; then
        log_info "Running in dry-run mode - no files will be downloaded"
    fi

    # Check authentication
    if [[ "$SKIP_AUTH_CHECK" != true ]]; then
        check_auth || {
            log_warning "Continuing without authentication - may fail for private repositories"
        }
    fi

    # Create directories
    create_directories

    # Download files
    local failed=0
    local success=0

    echo ""
    log_info "Starting downloads..."
    echo ""

    for key in "${!FILES[@]}"; do
        if download_file "$key"; then
            ((++success))
        else
            ((++failed))
        fi
    done

    # Summary
    echo ""
    echo "=========================================="
    echo "Download Summary"
    echo "=========================================="
    log_success "Successfully downloaded: $success files"
    if [[ $failed -gt 0 ]]; then
        log_error "Failed to download: $failed files"
        log_info "If authentication is required, set GITLAB_TOKEN and try again:"
        log_info "  export GITLAB_TOKEN=your-token-here"
        log_info "  ./download_test_data.sh"
        exit 1
    else
        log_success "All files downloaded successfully!"
    fi
}

# Run main function
main
