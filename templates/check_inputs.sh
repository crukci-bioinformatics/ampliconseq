#!/bin/bash

set -e -o pipefail

check_sample_sheet.R \
    --input !{sample_sheet} \
    --output !{checked_sample_sheet}

check_specific_variants.R \
    --input !{specific_variants} \
    --samples !{checked_sample_sheet} \
    --amplicons !{amplicons} \
    --output !{checked_specific_variants}

check_blacklisted_variants.R \
    --input !{blacklisted_variants} \
    --output !{checked_blacklisted_variants}

